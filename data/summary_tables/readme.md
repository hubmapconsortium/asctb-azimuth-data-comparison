### Azimuth and ASCT+B Data summary statistics

The corresponding summary file for each organ's celltype-counts is stored in this
directory as ```<organ_name>.celltype_stats.csv```.

The overall summary file for ASCTB master data, and Azimuth reference data ingested 
is stored in ```data/summary_tables``` directory as ```<organ_name>.celltype_stats.csv```. The ASCTB summary file for all organs
will contain the number-of-matches-found for CellTypes as well as Biomarkers, between ASCTB and Azimuth.

The R-package openxlsx has an issue, while writing HYPERTEXT formula to file. For now implemented hardcoded formula, so to use it:
1. Double click on the HYPERTEXT formula.

2. Click on another cell. This forces Excel to recognize the HYPERTEXT formula.

3. Click on hyperlink to go to a sheet.