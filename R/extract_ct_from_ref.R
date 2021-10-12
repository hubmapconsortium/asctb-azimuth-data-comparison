# DOCSTRING:      Updated the script on 09/17/2021 for the purposes of issue #1: Aggregated counts mismatch
#                 09/21/2021 Enhancement of issue #2: Pipeline needs to generate a summary of stats for Azimuth references
#                 09/23/2021 Enhancement of issue #2: Generate stats for counts of each celltype, within an organ. Still awaiting confirmation for computation-logic.
#                 10/12/2021 New pipeline integrated to ingest Google-Sheets of ASCT+B V1.1. master datasets and generate their summaries.
#                             Refined the logic for count-generation in most places.
#                             Also modularized the code into utility_functions.R so that this script only has functions relevant to summary-generation and high-level pipeline functionality.
# AUTHOR:         Darshal Shetty/ Vikrant Deshpande/ Amber Ramesh

library(Seurat)
library(rjson)
library(httr)
library(gsheet)

source('R/utility_extract_ct_from_json.R')
source('R/utility_functions.R')
source('R/summary_computation_functions.R')



# Initialization for [Azimuth reference vs ASCTB master] stats generation
azimuth_organ_stats_cols <- c("Organ", "Num.Unique.Cell.Types", "Num.Unique.CT.Ontology.IDs", "Num.Total.Cells", "Num.Annotation.Levels", "Num.Unique.Biomarkers")
azimuth_organ_stats <- create_new_df(azimuth_organ_stats_cols)
azimuth.entire_set_of_biomarkers <- NA

asctb_organ_stats_cols <- c("Organ", "Num.Unique.Anatomical.Structs", "Num.Unique.Cell.Types", "Num.Unique.CT.Ontology.IDs", "Num.Matching.CT.Ontology.IDs", "Num.Unique.Biomarkers", "Num.Matching.Biomarkers")
asctb_organ_stats <- create_new_df(asctb_organ_stats_cols)
asctb.entire_set_of_biomarkers <- NA

AZIMUTH.REFERENCE_RDS_DIR <- "data/azimuth_references/"
ASCTB_TARGET_DIR <- "data/asctb_formatted_azimuth_data/"
SUMMARIES_DIR <- "data/summary_tables/"
AZIMUTH.ANNOTATION_FILES_BASE_URL <- 'https://raw.githubusercontent.com/satijalab/azimuth_website/master/static/csv/'
CONFIGS <- rjson::fromJSON(file = 'data/azimuth_asctb_comparison_config.json')$references


for (config in CONFIGS) {
  
  cat(paste0("\n\nInitiating the ingestion for ",config$name," Azimuth reference..."))
  azimuth.entire_set_of_biomarkers <- c()
  asctb.entire_set_of_biomarkers <- c()
  
  # Create an ASCT-B format table from this organ's Azimuth reference
  asct_table <- switch(
    config$mode %||% '',
    'nested-json' = { extract_ct_from_json(config$url) },
    process_config_for_azimuth(config))
  
  # Wrangle the Azimuth dataset to derive summary stats
  process_azimuth_ref_dataset_summary(config, asct_table, azimuth_organ_stats)
  
  
  cat(paste0("\nInitiating the ingestion for ",config$name," ASCT+B master-tables..."))
  # Pull the Master table from this organ's Google-Sheet
  asctb_master_table <- get_asctb_master_table_content(config)
  
  # Wrangle the ASCT+B dataset to derive summary stats, or just add a dummy entry when no ASCTB-Master table
  suppressWarnings(
      msg <- process_asctb_master_dataset_summary(config=config, asctb_master_table=asctb_master_table, asct_table_derived_from_azimuth=asct_table)
    ,   classes="warning")
  cat(msg)
  
  # Finally, write the Azimuth dataset formatted as per the ASCTB structure for usability on CCF-reporter
  suppressWarnings(
    write_asctb_structure(config$name, asct_table)
    , classes="warning")
}

azimuth_organ_stats <- azimuth_organ_stats[order(azimuth_organ_stats$Organ),]
asctb_organ_stats <- asctb_organ_stats[order(asctb_organ_stats$Organ),]


# finally write the All-Organs ASCTB vs Azimuth stats into CSV files.
write_df_to_csv(azimuth_organ_stats, paste0(SUMMARIES_DIR,"Azimuth.All_organs.stats.csv"))
write_df_to_csv(asctb_organ_stats, paste0(SUMMARIES_DIR,"ASCTB.All_organs.stats.csv"))
