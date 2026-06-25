import pandas as pd
import brightway2 as bw
import logging
import uuid

class Metalinvent:

    def __init__(self,bw_project,ei_db_name,method,new_bio_name,metalinvent_db_name,path_file,df_cf_iwp):

        """

        Args:
            bw_project: Name of your brightway project
            ei_db_name: Name of ecoinvent cutoff 3.12 database in your bw project
            method: "Method_1" or "Method_2"
            new_bio_name: User defined name of new biosphere database to store missing elementary flows in biosphere3
            metalinvent_db_name: User defined name of adjusted version of ecoinvent cutoff 3.12
            path_file: Path to Excel file provided by author of this repo determining missing elementary flows in some activities
            df_cf_iwp: Excel file with ACP and RESEDA characterization factors of dissipative (within IMPACT World+ LCIA method)
        """
        self.bw_project = bw_project
        self.ei_db_name = ei_db_name
        self.method = method
        if self.method not in ["Method_1","Method_2"]:
            raise ValueError("Method must be either 'Method_1' or 'Method_2'")
        self.new_bio_name = new_bio_name
        self.metalinvent_db_name = metalinvent_db_name
        self.ei_adj_dict = {}
        self.biosphere_resources_dict = {}
        self.path_file = path_file
        df_change = pd.read_excel(self.path_file,sheet_name="missing_amounts").fillna(0)
        self.BtHav = pd.read_excel(self.path_file,sheet_name="BtH_av",index_col=0).fillna(0)
        self.deposit_types = [x for x in list(set(self.BtHav.Deposit)) if x!=0]
        self.elements_names = pd.read_excel(self.path_file,sheet_name="elements_names")
        self.df_change = df_change[df_change.Method == self.method]
        self.df_cf_iwp = df_cf_iwp
        self.column_name = {"Adaptation to resources services loss (beta)":"ACP CF value",
                            "Resources services deficit (beta)":"RESEDA CF value"}
        self.impact_cat = ["Adaptation to resources services loss (beta)","Resources services deficit (beta)"]
        self.comment_extraction = "; Metalinvent adaptations: Missing extraction flows added as per "+self.method+" in metalinvent tool."
        self.comment_dissipation = " Missing dissipative flows added as per "+self.method+" in metalinvent tool."
        self.df_missing_flows = pd.DataFrame(
            columns=["Elem flow name", "Compartment_iw", "Sub-compartment_iw", "code"])

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

        self.mineral_host = {"spodumene": {"Spodumene":{"Amount":0.0373,"Short":"Li"}},
                             "magnesium sulfate": {"Kieserite":{"Amount":0.1746,"Short":"Mg"}}}
        ## Li content in Spodumene from Table 1 of https://doi.org/10.3390/min9060334
        ## Mg content in Kieserite from https://seprm.com/wp-content/uploads/2020/11/Bimpilas-G.-M.-G.N.-Anastassakis.-Magnesite-beneficiation-methods-a-review.pdf



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

        self.subcomp_lookup = (
            self.bio3_flows.set_index("code")["Subcompartment"]
            .to_dict()
        )

        self.market_price_lookup = (
            self.elements_names
            .set_index("Long_Name")["market_price"]
            .to_dict()
        )

    def launch_operations(self):
        if self.method == "Method_2":
            self.logger.info("Building df change compliant with Method 2...")
            self.build_df_change_method2()
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

    def determine_host(self,process_code,ei_dict):
        print(process_code)
        nat_resource_codes = [
            exc["flow"]
            for exc in ei_dict[(self.ei_db_name, process_code)]["exchanges"]
            if self.subcomp_lookup.get(exc["flow"]) == "in ground"
        ]
        nat_resource_codes = list(set(nat_resource_codes))
        pass_nat_res = list(self.bio3_flows[self.bio3_flows.loc[:,"Elem flow name"].isin(list(self.elements_names.Long_Name))].code)
        nat_resource_codes = [x for x in nat_resource_codes if x in pass_nat_res]
        contributions = [
            (
                exc["name"],
                exc["amount"] * self.market_price_lookup.get(exc["name"], 0)
            )
            for exc in ei_dict[(self.ei_db_name, process_code)]["exchanges"]
            if exc["flow"] in nat_resource_codes
        ]
        if not contributions:
            return 0, 0

        total_contrib = sum(contrib for _, contrib in contributions)

        max_name, max_contrib = max(contributions, key=lambda x: x[1])
        host_short = self.elements_names[self.elements_names.Long_Name == max_name].Short_Name.iloc[0]
        share = max_contrib / total_contrib
        return host_short,share



    def build_df_change_method2(self):
        ei_dict = self.ei_db.load()
        df_change_mining = self.df_change[self.df_change.Analysis=="Mining"]

        for index in df_change_mining.index:
            process_code = self.ei_flows_with_codes[(self.ei_flows_with_codes.loc[:,"Product"]==self.df_change.loc[index,"reference product"])&
                                                    (self.ei_flows_with_codes.loc[:,"Activity"]==self.df_change.loc[index,"Activity"])&
                                                    (self.ei_flows_with_codes.loc[:,"Location"]==self.df_change.loc[index,"Location"])].code.iloc[0]
            if self.df_change.loc[index,"reference product"] not in self.mineral_host.keys():
                host_short = self.determine_host(process_code,ei_dict)[0]
                if not pd.isna(host_short) and host_short != 0:
                    host_long = self.elements_names[self.elements_names.Short_Name==host_short].Long_Name.iloc[0]
                    code = self.bio3_flows[(self.bio3_flows.loc[:,"Elem flow name"]==host_long)&
                                           (self.bio3_flows.loc[:,"Compartment"]=="natural resource")&
                                            (self.bio3_flows.loc[:, "Subcompartment"] == "in ground")].code.iloc[0]
                    amount_host = [exc["amount"] for exc in ei_dict[(self.ei_db_name,process_code)]["exchanges"] if exc["flow"]==code][0]
                    if len(self.BtHav.loc[host_short])>1:
                        if self.df_change.loc[index,"reference product"] in self.deposit_types:
                            deposit = self.df_change.loc[index,"reference product"]
                        else:
                            deposit = 0
                    else:
                        deposit=0
                    depo_row = self.BtHav[(self.BtHav.index == host_short) & (self.BtHav["Deposit"] == deposit)]
                    byproducts = [x for x in depo_row.columns[depo_row.loc[host_short] != 0].tolist() if len(x)<3]
                    for e in byproducts:
                        amount_byproduct = amount_host * depo_row.loc[host_short, e]
                        new_row = (
                            self.df_change.loc[index]
                            .copy()
                        )
                        new_row.update({
                            "Substance": e,
                            "Host": host_long,
                            "Substance_long_name":
                                self.elements_names[
                                    self.elements_names.Short_Name == e
                                    ].Long_Name.iloc[0],
                            "Missing extraction": amount_byproduct,
                            "Missing dissipation": amount_byproduct,
                            "Analysis":"Mining",
                        })

                        self.df_change.loc[len(self.df_change)] = new_row

            else:
                if self.df_change.loc[index, "reference product"] in self.deposit_types:
                    deposit = self.df_change.loc[index, "reference product"]
                else:
                    deposit = 0
                host = list(self.mineral_host[self.df_change.loc[index,"reference product"]].keys())[0]
                host_element = self.mineral_host[self.df_change.loc[index,"reference product"]][host]["Short"]
                depo_row = self.BtHav[(self.BtHav.index == host_element) & (self.BtHav["Deposit"] == deposit)]
                host_long = self.elements_names[self.elements_names.Short_Name == host_element].Long_Name.iloc[0]
                code = self.bio3_flows[(self.bio3_flows.loc[:, "Elem flow name"] == self.df_change.loc[index, "Host"]) &
                                       (self.bio3_flows.loc[:, "Compartment"] == "natural resource")&
                                        (self.bio3_flows.loc[:, "Subcompartment"] == "in ground")].code.iloc[0]
                amount_host_min = [exc["amount"] for exc in ei_dict[(self.ei_db_name, process_code)]["exchanges"] if
                               exc["flow"] == code][0]
                byproducts = [x for x in depo_row.columns[depo_row.loc[host_element] != 0].tolist() if len(x)<3]
                for e in byproducts:
                    amount_host = amount_host_min*self.mineral_host[self.df_change.loc[index,"reference product"]][host]["Amount"]
                    amount_byproduct = amount_host*depo_row.loc[host_element,e]
                    new_row = (
                        self.df_change.loc[index]
                        .copy()
                    )
                    new_row.update({
                        "Host": host_long,
                        "Substance":e,
                        "Substance_long_name":
                            self.elements_names[
                                self.elements_names.Short_Name == e
                                ].Long_Name.iloc[0],
                        "Missing extraction": amount_byproduct,
                        "Missing dissipation": amount_byproduct,
                        "Analysis": "Mining",
                    })

                    self.df_change.loc[len(self.df_change)] = new_row

    def complete_LCIA_methods(self):
        iw_methods = [method for method in bw.methods if "impact world+" in " ".join(method).lower()]
        for method in self.impact_cat:
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
                    print(((self.new_bio_name, self.df_new_bio_flows.loc[i, "code"]), self.df_new_bio_flows.loc[i, self.column_name[method]]))
                    try:
                        if self.df_new_bio_flows.loc[i, self.column_name[method]] !=0:
                            dict_method.append(((self.new_bio_name, self.df_new_bio_flows.loc[i, "code"]),
                                                self.df_new_bio_flows.loc[i, self.column_name[method]]))
                    except IndexError:
                        pass
                new_method.write(dict_method)

    def find_missing_flows(self):
        count_miss_flow = 0
        for i in self.df_change.index:
            qt_miss_diss = self.df_change.loc[i, "Missing dissipation"]
            if qt_miss_diss > 0:
                row_missing_flow = {"Elem flow name": self.df_change.loc[
                                                          i, "Substance_long_name"] + ", dissipative flow, to the environment",
                                    "Compartment_iw": "unspecified",
                                    "Sub-compartment_iw": "unspecified","code":uuid.uuid4().hex}
                self.df_missing_flows = pd.concat(
                    [self.df_missing_flows, pd.DataFrame(row_missing_flow, index=[count_miss_flow])])
        self.df_missing_flows = self.df_missing_flows.drop_duplicates(subset=["Elem flow name","Compartment_iw","Sub-compartment_iw"])
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
            self.biosphere_resources_dict[(self.new_bio_name, self.df_missing_flows.loc[i, "code"])] = {
                "name": self.df_missing_flows.loc[i, "Elem flow name"],
                "unit": "kilogram",
                "type": "biosphere",
                "categories": (self.df_missing_flows.loc[i, "Compartment_iw"],),
                "code": self.df_missing_flows.loc[i, "code"]
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
            for method in self.impact_cat:
                if "Aluminium" in self.df_new_bio_flows.loc[i, "Elem flow name"]:
                    if len(self.df_cf_iwp[(self.df_cf_iwp.loc[:, "Impact category"] == method) & (
                            self.df_cf_iwp.loc[:, "Elem flow name"] == "Aluminum, dissipative flow, to the environment")].loc[:,"CF value"]) > 0:
                        self.df_new_bio_flows.loc[i, self.column_name[method]] = \
                        self.df_cf_iwp[(self.df_cf_iwp.loc[:, "Impact category"]==method)&(self.df_cf_iwp.loc[:, "Elem flow name"] == "Aluminum, dissipative flow, to the environment")].loc[:,
                        "CF value"].iloc[0]
                    else:
                        self.df_new_bio_flows.loc[i, self.column_name[method]] = 0
                elif "Average plastic" in self.df_new_bio_flows.loc[i, "Elem flow name"]:
                    if method == self.impact_cat[0]:
                        self.df_new_bio_flows.loc[i, self.column_name[method]] = 3.15
                    else:
                        self.df_new_bio_flows.loc[i, self.column_name[method]] = 0
                else:
                    if len(self.df_cf_iwp[(self.df_cf_iwp.loc[:, "Impact category"]==method)&(self.df_cf_iwp.loc[:, "Elem flow name"] == self.df_new_bio_flows.loc[i, "Elem flow name"])].loc[:,
                    "CF value"])>0:
                        self.df_new_bio_flows.loc[i, self.column_name[method]] = \
                        self.df_cf_iwp[(self.df_cf_iwp.loc[:, "Impact category"]==method)&(self.df_cf_iwp.loc[:, "Elem flow name"] == self.df_new_bio_flows.loc[i, "Elem flow name"])].loc[:,
                        "CF value"].iloc[0]
                    else:
                        self.df_new_bio_flows.loc[i, self.column_name[method]] = 0

        self.df_new_bio_flows = self.df_new_bio_flows.drop_duplicates()
        self.df_new_bio_flows = self.df_new_bio_flows.reset_index().drop(columns=["index"])

    def copy_ei_db(self):
        """
        Here we copy ecoinvent and adjust the name of db_name in all exchanges between ecoinvent nodes.
        Then, we add missing extraction flows and missing dissipative flows
        Returns:

        """
        if self.metalinvent_db_name in list(bw.databases):
            del bw.databases[self.metalinvent_db_name]

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
            location = self.df_change.loc[i, "Location"]
            name = self.df_change.loc[i, "Activity"]
            refProduct = self.df_change.loc[i, "reference product"]
            analysis = self.df_change.loc[i, "Analysis"]
            if analysis == "EOL":
                self.process_codes_list = list(self.ei_flows_with_codes[
                                              (self.ei_flows_with_codes.loc[:, "Product"] == refProduct) & (
                                                          self.ei_flows_with_codes.loc[:, "Activity"] == name)].loc[:,
                                          "code"])
            elif analysis == "Mining":
                self.process_codes_list = list(self.ei_flows_with_codes[
                                              (self.ei_flows_with_codes.loc[:, "Location"] == location) & (
                                                          self.ei_flows_with_codes.loc[:, "Product"] == refProduct) & (
                                                          self.ei_flows_with_codes.loc[:, "Activity"] == name)].loc[:,
                                          "code"])
            if self.method == "Method_2":
                for code in self.process_codes_list:
                    for exc in self.ei_adj_dict[(self.metalinvent_db_name, code)]['exchanges']:
                        if "tailing" in exc["name"]:
                            exc["amount"] = 0
                            exc["comment"] += " Set to zero as per Method 2 of Metalinvent tool. Metal emissions from tailings are added as direct emissions in mine operation activity"
            elem_flow_name = self.df_change.loc[i, "Substance_long_name"]
            qt_missing_ext = self.df_change.loc[i, "Missing extraction"]
            if qt_missing_ext > 0:
                compartment = "natural resource"
                for code in self.process_codes_list:
                    code_flow_list = self.bio3_flows[
                                    (self.bio3_flows.loc[:, "Elem flow name"] == elem_flow_name) & (
                                                self.bio3_flows.loc[:, "Compartment"] == compartment)].loc[:,
                                "code"]
                    if len(code_flow_list)==1:
                        code_flow = code_flow_list.iloc[0]
                        biosphere_db = "biosphere3"
                        self.ei_adj_dict[(self.metalinvent_db_name, code)]['exchanges'].append({
                            "flow": code_flow,
                            "type": "biosphere",
                            "name":elem_flow_name,
                            "amount": self.df_change.loc[i, "Missing extraction"],
                            "input": (biosphere_db, code_flow),
                            "output": (self.metalinvent_db_name, code),
                            "comment": f"Missing extraction flow added as per {self.method} in metalinvent tool"
                        })
                        if 'comment' in list(self.ei_adj_dict[(self.metalinvent_db_name, code)].keys()):
                            if self.comment_extraction not in self.ei_adj_dict[(self.metalinvent_db_name, code)]['comment']:
                                self.ei_adj_dict[(self.metalinvent_db_name, code)]['comment'] += self.comment_extraction
                        else:
                            self.ei_adj_dict[(self.metalinvent_db_name, code)]['comment'] = self.comment_extraction

            qt_missing_diss = self.df_change.loc[i, "Missing dissipation"]
            if qt_missing_diss > 0:
                compartment = "unspecified"
                elem_flow_name = self.df_change.loc[i, "Substance_long_name"] + ", dissipative flow, to the environment"
                for code in self.process_codes_list:
                    if elem_flow_name in list(self.df_new_bio_flows.loc[:, "Elem flow name"]):
                        biosphere_db = self.new_bio_name
                        code_flow_list = \
                        self.df_new_bio_flows[self.df_new_bio_flows.loc[:, "Elem flow name"] == elem_flow_name].loc[:, "code"]
                        print("code_flow = ", code_flow_list)
                    else:
                        biosphere_db = "biosphere3"
                        code_flow_list = self.bio3_flows[
                                        (self.bio3_flows.loc[:, "Elem flow name"] == elem_flow_name) & (
                                                    self.bio3_flows.loc[:, "Compartment"] == compartment)].loc[:,
                                    "code"]
                    if len(code_flow_list) == 1:
                        code_flow = code_flow_list.iloc[0]
                        self.ei_adj_dict[(self.metalinvent_db_name, code)]['exchanges'].append({
                            "flow": code_flow,
                            "type": "biosphere",
                            "name": elem_flow_name,
                            "amount": self.df_change.loc[i, "Missing dissipation"],
                            "input": (biosphere_db, code_flow),
                            "output": (self.metalinvent_db_name, code),
                            "comment": f"Missing dissipative flow added as per {self.method} in metalinvent tool"
                        })
                        if 'comment' in list(self.ei_adj_dict[(self.metalinvent_db_name, code)].keys()):
                            if self.comment_dissipation not in self.ei_adj_dict[(self.metalinvent_db_name, code)]['comment']:
                                self.ei_adj_dict[(self.metalinvent_db_name, code)]['comment'] += self.comment_dissipation
                        else:
                            self.ei_adj_dict[(self.metalinvent_db_name, code)]['comment'] = self.comment_dissipation



        ### STEP 3 : Writing new metalinvent db into project
        self.logger.info("Writing metalinvent db to bw...")
        bw.Database(self.metalinvent_db_name).write(self.ei_adj_dict)


