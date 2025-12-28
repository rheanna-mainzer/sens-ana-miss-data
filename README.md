# sens-ana-miss-data
R code to carry out the analysis reported in the following paper:

Sensitivity analysis for missing data in observational studies: A practical guide for planning, conducting and reporting.

MASTER.R will run the following in order:
- descr_LSAC.R - clean and describe LSAC data
- cr_Figure2.R - create Figure 2 in paper
- sens_params_cont.R - inform choice of sensitivity paramneter for analysis of continuous outcome
- an_LSAC_cont.R - conduct analysis for the continuous outcome
- sens_params_bin.R - inform choic of sensitivity parameters for analysis of binary outcome
- an_LSAC_bin.R - conduct analysis for the binary outcome
- an_LSAC_cont_int - conduct extended analysis of continuous outcome reported in Supplement

*The data that support the findings of this study are available from the Australian Government Department of Social Services (DSS). Restrictions apply to the availability of these data, which were used under license for this study. Data are available at https://dataverse.ada.edu.au with the permission of the DSS.
