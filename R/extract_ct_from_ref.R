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
    print(paste('Something went wrong while processing the config for:',config$name))
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
    print(paste('Something went wrong while processing the config for:',config$name))
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
                  entire_set_of_biomarkers <<-  c(entire_set_of_biomarkers,unlist(strsplit(ont_data$Markers,",")))

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
      print(paste('Something went wrong while processing the config for:',config$name))
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




process_azimuth_ref_dataset_summary <- function(config, asct_table, organ_stats){
  tryCatch({
    # C1: Organ Name
    organ.1 <- config$asctb_name
    
    # C2: Get the union of all "CT" columns in the ASCTB organ.csv file
    entire_set_of_cell_types <- asct_table[!(grepl("ID",colnames(asct_table)) | grepl("COUNT",colnames(asct_table)) | grepl("LABEL",colnames(asct_table)))]
    n_unique_cell_types.2 <- length(unique(as.vector(as.matrix(entire_set_of_cell_types))))
    
    # C3: Get the union of all "ID" columns in the ASCTB organ.csv file
    entire_set_of_ct_ontology_ids <- asct_table[grepl("ID",colnames(asct_table))]
    n_unique_ct_ontology_ids.3 <- length(unique(as.vector(as.matrix(entire_set_of_ct_ontology_ids))))
    
    # C4: Get the sum of all counts that were already generated using the groupby formula used in process_config()
    n_total_cell_types.4 <- sum(asct_table[grep('COUNT',colnames(asct_table))])
    
    # C5: Get the number of annotations for an organ directly from the json config file
    n_annotation_levels.5 <- length(config$new_cell_type_files)
    
    # C6: Get the global variable containing all biomarkers for an organ, that we set during file-processing
    n_unique_biomarkers.6 <- length(unique(entire_set_of_biomarkers))
    
    # Append (C1, C2, ... C6) to the existing organ_stats global var
    organ_stats <<- rbind(organ_stats, c( organ.1, n_unique_cell_types.2, n_unique_ct_ontology_ids.3, 
                      n_total_cell_types.4, n_annotation_levels.5, n_unique_biomarkers.6))
    colnames(organ_stats) <<- organ_stats_cols
  },
  error = function(e){
      print(paste('Something went wrong while generating the organ-level summary stats for:',config$name))
      print(e)
  })
}









# Main initialization of global data-structures to capture summaries of each organ
organ_stats_cols <- c("Organ", "Num.Unique.Cell.Types", "Num.Unique.CT.Ontology.IDs", "Num.Total.Cells", "Num.Annotation.Levels", "Num.Unique.Biomarkers")
organ_stats <- data.frame(matrix(ncol=length(organ_stats_cols), nrow=0, dimnames=list(NULL,organ_stats_cols)))
entire_set_of_biomarkers <- c()
REFERENCE_RDS_DIR <- "data/azimuth_references/"
SUMMARIES_DIR <- "data/azimuth_summary_tables/"
ASCTB_TARGET_DIR <- "data/asctb_tables/"
ANNOTATION_FILES_BASE_URL <- 'https://raw.githubusercontent.com/satijalab/azimuth_website/master/static/csv/'



# main loop runs for each organ in the JSON config
for (config in rjson::fromJSON(file = 'data/organ_data.json')$references) {
  print(config$name)
  entire_set_of_biomarkers <- c()
  asct_table <- switch(
    config$mode %||% '',
    'nested-json' = { extract_ct_from_json(config$url) },
    process_config(config)
  )
  process_azimuth_ref_dataset_summary(config, asct_table, organ_stats)
  write_asctb(config$name, asct_table)
}

# finally write the Organ-level summaries into a CSV file.
write.table(organ_stats,
            file = paste0(SUMMARIES_DIR,"overall_organs.stats.csv"),
            sep = ',',
            na = "",
            append = FALSE,
            row.names = FALSE,
            col.names = TRUE)