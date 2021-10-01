# UPDATES:        Updated the script on 09/17/2021 for the purposes of issue #1: Aggregated counts mismatch
#                 Updated on 09/21/2021 for enhancement of issue #2: Pipeline needs to generate a summary of stats for Azimuth references
#                 Updated on 09/23/2021 for enhancement of issue #2: Generate stats for counts of each celltype, within an organ
# AUTHOR:         Darshal Shetty/ Vikrant Deshpande/ Amber Ramesh

# PBMC :          Azimuth-backend file doesn't contain 'AS/1'='other', which exists in the Azimuth reference dataset
# Motor Cortex :  Azimuth-backend file doesn't contain 'AS/4'='exclude', which exists in the Azimuth reference dataset
library(Seurat)
library(rjson)
library(httr)
source('R/extract_ct_from_json.R')


generate_celltype_vs_count_table <- function(cell_types, body_organ){
  tryCatch({
    
    # todo: Correct this logic for each combination of cell-type across levels.
    # Convert the available CellType columns in the RDS reference into a flattened vector
    all_celltype_values <- as.vector(as.matrix(cell_types))
    
    # Generate the pivot-table which shows celltype, count(*) for the organ
    celltype_vs_counts <- as.data.frame(table(all_celltype_values))
    
    # Basic reformatting for standardization
    renamed_cols <- c('Cell.Types', 'Num.Cells')
    colnames(celltype_vs_counts) <- renamed_cols
    celltype_vs_counts['Organ'] <- body_organ
    celltype_vs_counts <- celltype_vs_counts[c('Organ',renamed_cols)]
    
    write.csv(celltype_vs_counts, paste0(SUMMARIES_DIR, body_organ,".celltype_stats.csv"), row.names = FALSE)
  },
  error = function(e){
    print(paste('Something went wrong while generating the Celltype-vs-Counts summary for:',config$name))
    print(e)
  })
}




process_reference <- function(organ_config) {
  tryCatch({
    body_organ <- organ_config$name
    cell_hierarchy_cols <- organ_config$cell_type_columns # Still require this field for parsing the RDS data
    ref_url <- organ_config$url # Azimuth reference RDS-file download URL
  
    # load input reference files locally if present, else download them from URL in config
    file_name <- paste0(REFERENCE_RDS_DIR, body_organ, ".Rds")
    if (!file.exists(file_name)) {
      GET(ref_url, write_disk(file_name))
    }
    ref <- readRDS(file_name)
  
    # get columns containing cell type data
    cell_types <- ref@meta.data[cell_hierarchy_cols];
    # rename columns according to ASCT+B table format
    final_column_names <- sapply(1:length(cell_hierarchy_cols),
                                 function (x)
                                   return(paste0("AS/",x)))
    names(cell_types) <- final_column_names
  
    # generate a separate summary table: celltype, count(*) for this organ
    generate_celltype_vs_count_table(cell_types, body_organ)
    
    
    # get unique rows along with their counts
    count_col <- paste0("AS/", length(cell_hierarchy_cols), "/COUNT")
    cell_types[count_col] <- 1;
    group_by_formula <- as.formula(
      paste(
        paste0('`',count_col,'` ~ '),
        paste(sapply(final_column_names,
                     function(col)
                       paste0("`",col,"`")),
              collapse = " + ")))
    cell_types <- aggregate( group_by_formula, data = cell_types, FUN = sum )
  
    return(cell_types)
  },
  error = function(e){
    print(paste('Something went wrong while processing the Azimuth reference for:',config$name))
    print(e)
  })
}




process_cell_type_data <- function(body_organ, cell_hierarchy_cols) {
  tryCatch({
    # Pull the cell-annotation data from the official Azimuth Website backend-repository
    return(sapply(1:length(cell_hierarchy_cols),
                function (i, levels) {
                  # append the individual csv file-name for the organ into the base-url of raw-github-csv files
                  ct_file <- paste0(ANNOTATION_FILES_BASE_URL, levels[i])
                  
                  # read file containing cell ontology data for each cell type column in a body organ reference file
                  ont_data <- read.csv(ct_file)

                  # extract ontology label and id from "OBO.Ontology.ID" column
                  df <- data.frame(
                    name = ont_data$Label, # column to join with reference dataframe column values
                    label = if("OBO.Ontology.ID" %in% names(ont_data))
                              gsub("^\\[(.*?)\\].*",
                                   "\\1",
                                   ont_data$OBO.Ontology.ID)
                            else rep(NA, nrow(ont_data)),
                    id = if("OBO.Ontology.ID" %in% names(ont_data))
                          gsub(".*http:\\/\\/purl\\.obolibrary\\.org\\/obo\\/(.*?)_(.*?)\\)$",
                               "\\1:\\2",
                               ont_data$OBO.Ontology.ID)
                          else rep(NA, nrow(ont_data)))
                  
                  # set global variable of all biomarkers, so that we can create the summary of this organ later
                  azimuth.entire_set_of_biomarkers <<-  c(azimuth.entire_set_of_biomarkers, unlist(strsplit(ont_data$Markers,",")))

                  # rename columns according to ASCT+B table format
                  names(df) <- c(paste0("AS/",i),
                                 paste0("AS/",i,"/LABEL"),
                                 paste0("AS/",i,"/ID"))
                  return(df)
                  },
                cell_hierarchy_cols,
                simplify = FALSE))
  },
  error = function(e){
    print(paste('Something went wrong while processing the cell-type annotation file for:',config$name))
    print(e)
  })
}




write_asctb <- function(body_organ, asctb_table) {
  tryCatch({
    # Simply adds a metadata section into the final file, and then appends the actual contents
    column_count <- ncol(asctb_table)
    header_rows <- matrix(c(paste(body_organ,
                                  "cell types as anatomical structures",
                                  "from Azimuth reference data"), rep(NA, column_count-1),
                            rep(NA, column_count),
                            "Author Name(s):", "Azimuth & MC-IU", rep(NA, column_count-2),
                            "Author ORCID(s):", rep(NA, column_count-1),
                            "Reviewer(s):", rep(NA, column_count-1),
                            "General Publication(s):", rep(NA, column_count-1),
                            "Data DOI:", rep(NA, column_count-1),
                            "Date:", format(Sys.time(), "%x"), rep(NA, column_count-2),
                            "Version Number:", "v1.0", rep(NA, column_count-2),
                            rep(NA, column_count)),
                         10,
                         column_count,
                         byrow = TRUE)
    write.table(header_rows,
                file =  paste0(ASCTB_TARGET_DIR, body_organ, ".csv"),
                sep = ',',
                na = "",
                row.names = FALSE,
                col.names = FALSE)
    write.table(asctb_table,
                file = paste0(ASCTB_TARGET_DIR, body_organ, ".csv"),
                sep = ',',
                na = "",
                append = TRUE,
                row.names = FALSE,
                col.names = TRUE)
    },
    error = function(e){
      print(paste('Something went wrong while writing the ASCTB table for:',config$name))
      print(e)
    })
}




process_config <- function(config) {
  tryCatch({
    # extract azimuth reference data
    reference_table <- process_reference(config)
  
    #pull the files for cell-type ontology data from the backend Azimuth website Repo directly
    ct_ontology_tables <- process_cell_type_data(config$name, config$new_cell_type_files)
  
    # map ontology ID and LABELS to reference cell types
    merged_data <- Reduce(merge, ct_ontology_tables, reference_table)
    
    # reorder columns
    column_order <- c(
      sapply(1:length(config$new_cell_type_files),
             function (n)
               c(paste0("AS/",n),
                 paste0("AS/",n,"/LABEL"),
                 paste0("AS/",n,"/ID"))),
      paste0("AS/",length(config$new_cell_type_files),"/COUNT"))
    
    # generate final CSV file
    return(merged_data[,column_order])
  },
  error = function(e){
    print(paste('Something went wrong while processing the config for:',config$name))
    print(e)
  })
}




process_azimuth_ref_dataset_summary <- function(config, asct_table, azimuth_organ_stats){
  tryCatch({
    asctb.formatted_azimuth_columns <- colnames(asct_table)
    
    # C1: Organ Name
    organ.1 <- config$asctb_name
    
    # C2: Get the union of all "CT" columns in the ASCTB organ.csv file
    entire_set_of_cell_types <- asct_table[grepl("AS/[0-9]$",asctb.formatted_azimuth_columns)]
    cleaned_set_of_cell_types <- get_cleaned_values_from_df(entire_set_of_cell_types)
    n_unique_cell_types.2 <- length(cleaned_set_of_cell_types)
    
    # C3: Get the union of all "ID" columns in the ASCTB organ.csv file
    entire_set_of_ct_ontology_ids <- asct_table[grepl("AS",asctb.formatted_azimuth_columns) & grepl("ID",asctb.formatted_azimuth_columns)]
    cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(entire_set_of_ct_ontology_ids)
    n_unique_ct_ontology_ids.3 <- length(cleaned_set_of_ct_ontology_ids)
    
    # C4: Get the sum of all counts that were already generated using the groupby formula used in process_config()
    n_total_cell_types.4 <- sum(asct_table[grep('COUNT',asctb.formatted_azimuth_columns)])
    
    # C5: Get the number of annotations for an organ directly from the json config file
    n_annotation_levels.5 <- length(config$new_cell_type_files)
    
    # C6: Get the global variable containing all biomarkers for an organ, that we set during file-processing
    azimuth.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(azimuth.entire_set_of_biomarkers)
    n_unique_biomarkers.6 <- length(azimuth.cleaned_set_of_biomarkers)
    
    # Append (C1, C2, ... C6) to the existing azimuth_organ_stats global var
    azimuth_organ_stats <<- rbind(azimuth_organ_stats, c( organ.1, n_unique_cell_types.2, n_unique_ct_ontology_ids.3, 
                                                          n_total_cell_types.4, n_annotation_levels.5, n_unique_biomarkers.6))
    colnames(azimuth_organ_stats) <<- azimuth_organ_stats_cols
  },
  error = function(e){
    print(paste('Something went wrong while generating the Azimuth organ-level summary stats for:',config$name))
    print(e)
  })
}














##################### Attempting to read Google Sheet data for ASCTB Master Tables #####################

#install.packages("gsheet")
library(gsheet)


get_master_table_content <- function(config){
  # Read the google-sheet tab url and return a dataframe after replacing the metainformation in top 10 rows. If URL not available return NA.
  tryCatch({
    
    url <- config$asctb_master_url
    if (url=="NA"){
      return (NA)
    }
    asctb.master.data <- gsheet2tbl(url)
    
    # Remove out the top 10 meta info rows
    asctb.master.data <- asctb.master.data[10:nrow(asctb.master.data),]
    
    # Set colnames to the first available row
    colnames(asctb.master.data) <- asctb.master.data[1,]
    
    # Remove out the top row which was just the colnames
    asctb.master.data <- asctb.master.data[2:nrow(asctb.master.data),]
    return (asctb.master.data)
    
    },
    error=function(e){
      print(paste('Something went wrong while generating the ASCTB master table for:',config$name))
      print(e)
    })
}


get_cleaned_values_from_df <- function(df){
  tryCatch({
    if (is.null(df) || (is.data.frame(df) && (nrow(df)*ncol(df)==0)) || sum(is.na(df)==nrow(df))){
      return (c())
    }else{
      # Apply some basic functions to the entire df and return the unique values without NA
      cleaned_df_values <- na.omit(unique(trimws(as.vector(as.matrix(df)))))
      return (cleaned_df_values[!grepl("any", cleaned_df_values)])
    }
  },
  error=function(e){
    print('Something went wrong while retrieving the cleaned values of dataframe')
    print(e)
  })
}





process_asctb_master_dataset_summary <- function(config, asctb_master_table, asctb_organ_stats, asct_table_derived_from_azimuth){
  tryCatch({
    asctb.master_columns <- colnames(asctb_master_table)
    
    # C1: Organ Name
    organ.1 <- config$asctb_name
    
    
    # C2: Get the union of all "CT" columns in the ASCTB organ.csv file
    # Edge cases: PBMC needs 'CT' column-values added in as well.
    #             Spleen needs 'A S/6' column-values added in as well. Counts still don't match
    #if (config$name %in% c("pbmc","spleen")){
    #  asctb.entire_set_of_cell_types <- asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns) | grepl("A S/[0-9]$",asctb.master_columns)]
    #}else{
    #  asctb.entire_set_of_cell_types <- asctb_master_table[grepl("AS/[0-9]$",asctb.master_columns)]
    #}
    asctb.entire_set_of_cell_types <- asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns)]
    asctb.cleaned_set_of_cell_types <- get_cleaned_values_from_df(asctb.entire_set_of_cell_types)
    n_unique_cell_types.2 <- length(asctb.cleaned_set_of_cell_types)
    
    
    # C3: Get the union of all "ID" columns in the ASCTB organ.csv file. Kidney has extra 'CT' column that messes up the counts
    #if (config$name=='pbmc' || config$name=='spleen'){
    #  asctb.entire_set_of_ct_ontology_ids <- asctb_master_table[(grepl("AS",asctb.master_columns) | grepl("A S",asctb.master_columns) | grepl("CT",asctb.master_columns)) & grepl("ID",asctb.master_columns)]
    #}else{
    #  asctb.entire_set_of_ct_ontology_ids <- asctb_master_table[(grepl("AS",asctb.master_columns)) & grepl("ID",asctb.master_columns)]
    #}
    asctb.entire_set_of_ct_ontology_ids <- asctb_master_table[(grepl("CT",asctb.master_columns)) & grepl("ID",asctb.master_columns)]
    asctb.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(asctb.entire_set_of_ct_ontology_ids)
    n_unique_ct_ontology_ids.3 <- length(asctb.cleaned_set_of_ct_ontology_ids)
    
    
    
    # C4: Get the intersection of CT_Ontology_IDs in the ASCTB organ.csv file vs Azimuth organ.csv file
    azimuth_colnames <-  colnames(asct_table_derived_from_azimuth)
    azimuth.asctb.entire_set_of_ct_ontology_ids <- asct_table_derived_from_azimuth[grepl("AS",azimuth_colnames) & grepl("ID",azimuth_colnames)]
    azimuth.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(azimuth.asctb.entire_set_of_ct_ontology_ids)
    
    asctb.entire_set_of_ct_ontology_ids <- asctb_master_table[grepl("CT",asctb.master_columns) & grepl("ID",asctb.master_columns)]
    asctb.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(asctb.entire_set_of_ct_ontology_ids)
    n_matching_ct_ontology_ids.4 <- length(intersect(asctb.cleaned_set_of_ct_ontology_ids, azimuth.cleaned_set_of_ct_ontology_ids))
    
    
    # C5: Get the union of all "BG" columns in the ASCTB organ.csv file
    asctb.entire_set_of_biomarkers <<- asctb_master_table[(grepl("BG",asctb.master_columns)) & !(grepl("ID",asctb.master_columns) | grepl("COUNT",asctb.master_columns) | grepl("LABEL",asctb.master_columns))]
    asctb.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(asctb.entire_set_of_biomarkers)
    n_unique_biomarkers.5 <- length(asctb.cleaned_set_of_biomarkers)
    
    
    # C6: Get the intersection of Biomarkers in the ASCTB organ.csv file vs Azimuth organ.csv file
    azimuth.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(azimuth.entire_set_of_biomarkers)
    n_matching_biomarkers.6 <- length(intersect(azimuth.cleaned_set_of_biomarkers, asctb.cleaned_set_of_biomarkers))
    
    
    # Append (C1, C2, ... C6) to the existing asctb_organ_stats global var
    asctb_organ_stats <<- rbind(asctb_organ_stats, c( organ.1, n_unique_cell_types.2, n_unique_ct_ontology_ids.3, 
                                                      n_matching_ct_ontology_ids.4, n_unique_biomarkers.5, n_matching_biomarkers.6 ))
    colnames(asctb_organ_stats) <<- asctb_organ_stats_cols
  },
  error = function(e){
    print(paste('Something went wrong while generating the ASCTB organ-level summary stats for:',config$name))
    print(e)
  })
}

















# Main initialization of global data-structures to capture summaries of each organ
azimuth_organ_stats_cols <- c("Organ", "Num.Unique.Cell.Types", "Num.Unique.CT.Ontology.IDs", "Num.Total.Cells", "Num.Annotation.Levels", "Num.Unique.Biomarkers")
azimuth_organ_stats <- data.frame(matrix(ncol=length(azimuth_organ_stats_cols), nrow=0, dimnames=list(NULL,azimuth_organ_stats_cols)))
azimuth.entire_set_of_biomarkers <- NA

asctb_organ_stats_cols <- c("Organ", "Num.Unique.Cell.Types", "Num.Unique.CT.Ontology.IDs", "Num.Matching.CT.Ontology.IDs", "Num.Unique.Biomarkers", "Num.Matching.Biomarkers")
asctb_organ_stats <- data.frame(matrix(ncol=length(azimuth_organ_stats_cols), nrow=0, dimnames=list(NULL,azimuth_organ_stats_cols)))
asctb.entire_set_of_biomarkers <- NA

REFERENCE_RDS_DIR <- "data/azimuth_references/"
SUMMARIES_DIR <- "data/azimuth_summary_tables/"
ASCTB_TARGET_DIR <- "data/asctb_tables/"
ANNOTATION_FILES_BASE_URL <- 'https://raw.githubusercontent.com/satijalab/azimuth_website/master/static/csv/'



configs <- rjson::fromJSON(file = 'data/organ_data.json')$references
for (config in configs) {
  print(config$name)
  azimuth.entire_set_of_biomarkers <- c()
  asctb.entire_set_of_biomarkers <- c()
  
  # Create an ASCT-B format table from this organ's Azimuth reference
  asct_table <- switch(
    config$mode %||% '',
    'nested-json' = { extract_ct_from_json(config$url) },
    process_config(config))
  
  # Wrangle the Azimuth dataset to derive summary stats
  process_azimuth_ref_dataset_summary(config, asct_table, azimuth_organ_stats)
  
  
  # Pull the Master table from this organ's Google-Sheet
  asctb_master_table <- get_master_table_content(config)
  
  # Wrangle the ASCT+B dataset to derive summary stats, or just add a dummy entry when no ASCTB-Master table
  suppressWarnings(
      if (!is.na.data.frame(asctb_master_table)){
        process_asctb_master_dataset_summary(config=config, asctb_master_table=asctb_master_table, asctb_organ_stats=asctb_organ_stats, asct_table_derived_from_azimuth=asct_table)
      }else{
        asctb_organ_stats <- rbind(asctb_organ_stats, c(config$asctb_name, rep(0,length(colnames(asctb_organ_stats))-1)))
        colnames(asctb_organ_stats) <- asctb_organ_stats_cols
      }
    ,   classes="warning")
  
  # Finally, write the Azimuth dataset formatted as per the ASCTB structure so that they are usable on CCF-reporter
  suppressWarnings(
      write_asctb(config$name, asct_table) 
  , classes="warning")
}

asctb_organ_stats <- asctb_organ_stats[order(asctb_organ_stats$Organ),]


# finally write the All-Organs summaries into a CSV file.
write.table(azimuth_organ_stats,
            file = paste0(SUMMARIES_DIR,"all_organs.stats.csv"),
            sep = ',',
            na = "",
            append = FALSE,
            row.names = FALSE,
            col.names = TRUE)