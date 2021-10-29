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
library(openxlsx)

source('R/utility_extract_ct_from_json.R')
source('R/utility_functions.R')
source('R/summary_computation_functions.R')



# Initialization for [Azimuth reference vs ASCTB master] stats generation
azimuth_organ_stats_cols <- c("Organ", "AZ.Unique.Cell.Types", "AZ.Unique.CT.IDs", "AZ.Total.Cells", "AZ.Annotation.Levels", "AZ.Unique.Biomarkers", "Raw.Organ.Name")
azimuth_organ_stats <- create_new_df(azimuth_organ_stats_cols)
azimuth.entire_set_of_biomarkers <- NA

asctb_organ_stats_cols <- c("Organ", "ASCTB.Unique.Cell.Types", "ASCTB.Unique.CT.IDs", "Matching.CT.IDs", "ASCTB.Unique.Biomarker.Genes", "ASCTB.Unique.Biomarker.Prots", "Matching.Biomarkers")
asctb_organ_stats <- create_new_df(asctb_organ_stats_cols)
asctb.entire_set_of_biomarkers <- NA


AZIMUTH.REFERENCE_RDS_DIR <- "data/azimuth_references/"
ASCTB_TARGET_DIR <- "data/asctb_formatted_azimuth_data/"
SUMMARIES_DIR <- "data/summary_tables/"
STAGING_DIR <- "data/staging_area/"
AZIMUTH.ANNOTATION_FILES_BASE_URL <- 'https://raw.githubusercontent.com/satijalab/azimuth_website/master/static/csv/'
CONFIGS <- rjson::fromJSON(file = 'data/azimuth_asctb_comparison_config.json')$references
BIOMARKER_NAME_VS_ID_CACHE <- 'data/biomarker_name_vs_id_cached.csv'

BIOMARKER_NAME_VS_ID_MAPPING <- as.data.frame(read.csv(BIOMARKER_NAME_VS_ID_CACHE, na.string=c("NA", "NULL"), encoding="UTF-8"))
# config <- CONFIGS[[2]]

for (config in CONFIGS) {
  
  cat("\n\n\n\nInitiating the ingestion for ",config$name," Azimuth reference...")
  azimuth.entire_set_of_biomarkers <- c()
  asctb.entire_set_of_biomarkers <- c()
  
  # Create an ASCT-B format table from this organ's Azimuth reference
  asct_table <- switch(
    config$mode %||% '',
    'nested-json' = { extract_ct_from_json(config$url) },
    process_config_for_azimuth(config))
  
  # Wrangle the Azimuth dataset to derive summary stats
  process_azimuth_ref_dataset_summary(config, asct_table, azimuth_organ_stats)
  
  
  cat("\nInitiating the ingestion for ",config$name," ASCT+B master-tables...")
  # Pull the Master table from this organ's Google-Sheet
  asctb.file_path <- get_asctb_master_table_content(config)
  
  # Wrangle the ASCT+B dataset to derive summary stats, or just add a dummy entry when no ASCTB-Master table
  suppressWarnings(
      msg <- process_asctb_master_dataset_summary(config=config, file_path=asctb.file_path, asct_table_derived_from_azimuth=asct_table)
    ,   classes="warning")
  cat(msg)
  
  # Finally, write the Azimuth dataset formatted as per the ASCTB structure for usability on CCF-reporter
  suppressWarnings( write_asctb_structure(config$name, asct_table) , classes="warning")
}

# finally write the All-Organs ASCTB vs Azimuth stats into CSV files.
create_combined_summaries(asctb_organ_stats, azimuth_organ_stats)