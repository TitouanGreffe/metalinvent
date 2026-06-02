import pandas as pd
import brightway2 as bw
import logging
import uuid

class Metalinvent:

    def __init__(self,bw_project,ei_db_name,new_bio_name,metalinvent_db_name,df_change,df_cf_iwp):

        self.bw_project = bw_project
        self.ei_db_name = ei_db_name
        self.new_bio_name = new_bio_name
        self.metalinvent_db_name = metalinvent_db_name
        self.ei_adj_dict = {}
        self.biosphere_resources_dict = {}
        self.df_change = df_change
        self.df_cf_iwp = df_cf_iwp
        self.df_missing_flows = pd.DataFrame(
            columns=["Elem flow name", "Compartment_iw", "Sub-compartment_iw", "ACP CF value"])

        # set up logging tool
        self.logger = logging.getLogger('Metalinvent')
        self.logger.setLevel(logging.INFO)
        self.logger.handlers = []
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        self.logger.propagate = False




        bw.projects.set_current(self.bw_project)
        if "biosphere3" not in list(bw.databases):
            self.logger.info("biosphere3 is missing...")
        elif self.ei_db_name not in list(bw.databases):
            self.logger.info(str(self.ei_db_name)+" is missing...")
        else:
            self.ei_db = bw.Database(self.ei_db_name)
            self.bio = bw.Database("biosphere3")

        self.bio3_flows = pd.DataFrame(
            [
                (
                    i.as_dict()['name'],
                    i.as_dict()['categories'][0],
                    i.as_dict()['categories'][1] if len(i.as_dict()['categories']) == 2 else "",
                    i.as_dict()['code'],
                    i.as_dict()['unit']
                )
                for i in self.bio
            ],
            columns=['Elem flow name', 'Compartment', 'Subcompartment', 'code', 'unit']
        )
        self.ei_flows_with_codes = (
            pd.DataFrame(
                [(i.as_dict()['reference product'], i.as_dict()['name'], i.as_dict()['location'],
                  i.as_dict()['code'], i.as_dict()['unit'], self.get_classification(i.as_dict())[0],
                  self.get_classification(i.as_dict())[1])

                 for i in self.ei_db],
                columns=['Product', 'Activity', 'Location', 'code', 'unit', "ISIC sector", "CPC code"])
        )

    def launch_operations(self):
        self.logger.info("Finding missing flows...")
        self.find_missing_flows()
        self.logger.info("Creating new biosphere...")
        self.create_new_biosphere()
        self.logger.info("Writing new biosphere...")
        self.write_new_biosphere(self.biosphere_resources_dict)
        self.logger.info("Building df new biosphere ...")
        self.df_new_biosphere()
        self.logger.info("Copy ecoinvent, adjust and write in bw ...")
        self.copy_ei_db()
        self.logger.info("Loading LCIA methods and add CFs of newly added elementary flows..")
        self.complete_LCIA_methods()

    def complete_LCIA_methods(self):
        impact_cat = ["Adaptation to resources services loss (beta)"]
        iw_methods = [method for method in bw.methods if "impact world+" in " ".join(method).lower()]
        for method in impact_cat:
            method_bw = [m for m in iw_methods if "IMPACT World+ Midpoint 2.2.1 for ecoinvent v3.12" in m[0] if
                         method in m[2]]
            if len(method_bw)==0:
                self.logger.info(method+" is missing... Add the method into project")
            else:
                self.logger.info("Completing LCIA method with new elementary flows")
                new_method = bw.Method(method_bw[0])
                # register the new method
                new_method.register()
                # set its unit
                new_method.metadata["unit"] = bw.Method(method_bw[0]).metadata["unit"]
                dict_method = bw.Method(method_bw[0]).load()
                for i in self.df_new_bio_flows.index:
                    print(((self.new_bio_name, self.df_new_bio_flows[(self.df_new_bio_flows.loc[:, "Elem flow name"] ==
                                                                     self.df_new_bio_flows.loc[i, "Elem flow name"]) & (
                                                                                self.df_new_bio_flows.loc[:, "Compartment"] ==
                                                                                self.df_new_bio_flows.loc[i, "Compartment"])].loc[
                                                     :, "code"].iloc[0]), self.df_new_bio_flows.loc[i, 'ACP CF value']))
                    try:
                        dict_method.append(((self.new_bio_name, self.df_new_bio_flows[(self.df_new_bio_flows.loc[:,
                                                                                      "Elem flow name"] ==
                                                                                      self.df_new_bio_flows.loc[
                                                                                          i, "Elem flow name"]) & (
                                                                                                 self.df_new_bio_flows.loc[:,
                                                                                                 "Compartment"] ==
                                                                                                 self.df_new_bio_flows.loc[
                                                                                                     i, "Compartment"])].loc[
                                                                      :, "code"].iloc[0]),
                                            self.df_new_bio_flows.loc[i, 'ACP CF value']))
                    except IndexError:
                        pass
                new_method.write(dict_method)

    def find_missing_flows(self):
        count_miss_flow = 0
        for i in self.df_change.index:
            qt_miss_diss = self.df_change.loc[i, "Missing dissipation"]
            if qt_miss_diss > 0:
                row_missing_flow = {"Elem flow name": self.df_change.loc[
                                                          i, "Substance long name"] + ", dissipative flow, to the environment",
                                    "Compartment_iw": "unspecified",
                                    "Sub-compartment_iw": "unspecified", "ACP CF value": self.df_change.loc[i, "ACP CF"]}
                self.df_missing_flows = pd.concat(
                    [self.df_missing_flows, pd.DataFrame(row_missing_flow, index=[count_miss_flow])])
        self.df_missing_flows = self.df_missing_flows.drop_duplicates()
        self.df_missing_flows = self.df_missing_flows.reset_index().drop(columns=["index"])


    def get_classification(self,dict_process):
        if "classifications" not in dict_process:
            return ("", "")
        else:
            isic_sector = ""
            cpc_code = ""
            for c in dict_process["classifications"]:
                if c[0] == "ISIC rev.4 ecoinvent":
                    isic_sector = c[1]
                if c[0] == "CPC":
                    cpc_code = c[1]
            return (isic_sector, cpc_code)

    def create_new_biosphere(self):
        for i in self.df_missing_flows.index:
            code = uuid.uuid4().hex
            self.biosphere_resources_dict[(self.new_bio_name, code)] = {
                "name": self.df_missing_flows.loc[i, "Elem flow name"],
                "unit": "kilogram",
                "type": "biosphere",
                "categories": (self.df_missing_flows.loc[i, "Compartment_iw"],),
                "code": code
            }
    def write_new_biosphere(self,new_bio_dict):
        if self.new_bio_name in list(bw.databases):
            del bw.databases[self.new_bio_name]
        bw.Database(self.new_bio_name).write(new_bio_dict)

    def df_new_biosphere(self):
        self.df_new_bio_flows = (
            pd.DataFrame(
                [(self.biosphere_resources_dict[i]['name'], self.biosphere_resources_dict[i]['categories'][0],
                  self.biosphere_resources_dict[i]['categories'][1],
                  self.biosphere_resources_dict[i]['code'], self.biosphere_resources_dict[i]['unit'])
                 if len(self.biosphere_resources_dict[i]['categories']) == 2
                 else (self.biosphere_resources_dict[i]['name'], self.biosphere_resources_dict[i]['categories'][0], 'unspecified',
                       self.biosphere_resources_dict[i]['code'], self.biosphere_resources_dict[i]['unit'])
                 for i in self.biosphere_resources_dict.keys()],
                columns=['Elem flow name', 'Compartment', 'Sub-compartment', 'code', 'unit'])
        )
        for i in self.df_new_bio_flows.index:
            if "Aluminium" in self.df_new_bio_flows.loc[i, "Elem flow name"]:
                self.df_new_bio_flows.loc[i, "ACP CF value"] = \
                self.df_cf_iwp[self.df_cf_iwp.loc[:, "substances"] == "Aluminum, dissipative flow, to the environment"].loc[:,
                "ACP STEPS with C"].iloc[0]
            elif "Average plastic" in self.df_new_bio_flows.loc[i, "Elem flow name"]:
                self.df_new_bio_flows.loc[i, "ACP CF value"] = 3.15
            else:
                self.df_new_bio_flows.loc[i, "ACP CF value"] = \
                self.df_cf_iwp[self.df_cf_iwp.loc[:, "substances"] == self.df_new_bio_flows.loc[i, "Elem flow name"]].loc[:,
                "ACP STEPS with C"].iloc[0]
        self.df_new_bio_flows.drop_duplicates()
        self.df_new_bio_flows.reset_index().drop(columns=["index"])

    def copy_ei_db(self):
        """
        Here we copy ecoinvent and adjust the name of db_name in all exchanges between ecoinvent nodes.
        Then, we add missing extraction flows and missing dissipative flows
        Returns:

        """

        ## STEP 1: Loading ei_db and adjusting db key in new dictionary ei_adj_dict
        ei_dict = self.ei_db.load()
        self.ei_adj_dict = {}
        for (db, uid), value in ei_dict.items():
            if db == self.ei_db_name:
                db = self.metalinvent_db_name
            self.ei_adj_dict[(db, uid)] = value
            for exc in self.ei_adj_dict[(db, uid)]["exchanges"]:
                for key in ('input', 'output'):
                    if key in exc and exc[key][0] == self.ei_db_name:
                        exc[key] = (self.metalinvent_db_name, exc[key][1])

        ### STEP 2 : Adding missing extraction and dissipative flows
        for i in self.df_change.index:
            Location = self.df_change.loc[i, "Location"]
            name = self.df_change.loc[i, "Activity"]
            refProduct = self.df_change.loc[i, "reference product"]
            analysis = self.df_change.loc[i, "Analysis"]
            if analysis == "EOL":
                process_codes_list = list(self.ei_flows_with_codes[
                                              (self.ei_flows_with_codes.loc[:, "Product"] == refProduct) & (
                                                          self.ei_flows_with_codes.loc[:, "Activity"] == name)].loc[:,
                                          "code"])
            elif analysis == "Mining":
                process_codes_list = list(self.ei_flows_with_codes[
                                              (self.ei_flows_with_codes.loc[:, "Location"] == Location) & (
                                                          self.ei_flows_with_codes.loc[:, "Product"] == refProduct) & (
                                                          self.ei_flows_with_codes.loc[:, "Activity"] == name)].loc[:,
                                          "code"])
            elem_flow_name = self.df_change.loc[i, "Substance long name"]
            qt_missing_ext = self.df_change.loc[i, "Missing extraction"]
            if qt_missing_ext > 0:
                compartment = "natural resource"
                for code in process_codes_list:
                    code_flow = self.bio3_flows[
                                    (self.bio3_flows.loc[:, "Elem flow name"] == elem_flow_name) & (
                                                self.bio3_flows.loc[:, "Compartment"] == compartment)].loc[:,
                                "code"].iloc[0]
                    biosphere_db = "biosphere3"
                    self.ei_adj_dict[(self.metalinvent_db_name, code)]['exchanges'].append({
                        "flow": elem_flow_name,
                        "type": "biosphere",
                        "amount": self.df_change.loc[i, "Missing extraction"],
                        "input": (biosphere_db, code_flow),
                        "output": (self.metalinvent_db_name, code),
                        "comment": "Missing extraction flow added as per metalinvent tool"
                    })

            qt_missing_diss = self.df_change.loc[i, "Missing dissipation"]
            if qt_missing_diss > 0:
                compartment = "unspecified"
                elem_flow_name = self.df_change.loc[i, "Substance long name"] + ", dissipative flow, to the environment"
                for code in process_codes_list:
                    if elem_flow_name in list(self.df_new_bio_flows.loc[:, "Elem flow name"]):
                        biosphere_db = self.new_bio_name
                        code_flow = \
                        self.df_new_bio_flows[self.df_new_bio_flows.loc[:, "Elem flow name"] == elem_flow_name].loc[:, "code"].iloc[0]
                        print("code_flow = ", code_flow)
                    else:
                        biosphere_db = "biosphere3"
                        code_flow = self.bio3_flows[
                                        (self.bio3_flows.loc[:, "Elem flow name"] == elem_flow_name) & (
                                                    self.bio3_flows.loc[:, "Compartment"] == compartment)].loc[:,
                                    "code"].iloc[0]
                    self.ei_adj_dict[(self.metalinvent_db_name, code)]['exchanges'].append({
                        "flow": elem_flow_name,
                        "type": "biosphere",
                        "amount": self.df_change.loc[i, "Missing dissipation"],
                        "input": (biosphere_db, code_flow),
                        "output": (self.metalinvent_db_name, code),
                        "comment": "Missing dissipative flow added as per metalinvent tool"
                    })
        ### STEP 3 : Writing new metalinvent db into project
        if self.metalinvent_db_name in list(bw.databases):
            del bw.databases[self.metalinvent_db_name]
        else:
            bw.Database(self.metalinvent_db_name).write(self.ei_adj_dict)


