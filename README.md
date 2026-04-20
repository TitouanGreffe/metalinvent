# metalinvent
metalinvent: a brightway2 compliant tool to adjust abiotic resource flows in ecoinvent

This python tool allows to add missing extraction flows and dissipative flows in targeted ecoinvent datasets, namely ore mining activities and waste treatment activities.
When some dissipative flows are not included in the biosphere3 brightway database, it creates a new biosphere_resource database with new flows, include those flows in targeted activities and include IMPACT World+ version 2.2.1 characterization factors for flows in biosphere_resource.

# Methods

![metalinvent general method](doc/impact_pathway.png)

### Requirements
ecoinvent license to download ecoinvent 3.12 cutoff database.

# Author
Titouan Greffe (greffe.titouan@uqam.ca)
# Reference
Greffe, T., Agez M., Margni, M., Bulle, C. (2026) No abiotic resource flows, no impacts? Scrutinizing hidden potential impact of unreported extraction and dissipative flows in ecoinvent life cycle inventory database. _under preparation_
