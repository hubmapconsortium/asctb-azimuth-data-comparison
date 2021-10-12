# DOCSTRING:      Modularized the code into utility_functions.R so that extract_ct_from_ref.R only has logic for the entire pipeline.
#                 This script will have modular functions for 2 purposes:
#                  1. Ingest Azimuth reference data and convert into ASCTB format as follows:
#                      process_config_for_azimuth()-
#                        process_azimuth_reference() ingests the RDS Zenodo file.
#                        process_azimuth_annotation_celltype_data() ingests the annotation level-files for Azimuth.
#                        Merge both of these files to retrieve the Cell-type IDs from annotation level-files.
#                        write_asctb_structure() will write the merged dataframe with meta-information so that CCF-portal can visualize this structure.
#                  2. Ingest ASCTB Master-data from the Google-sheets available on CCF-portal as follows:
#                      get_asctb_master_table_content()-
#                        Reads the Google-sheets and removes the meta-information rows.
#
# AUTHOR:         Darshal Shetty/ Vikrant Deshpande/ Amber Ramesh


create_new_df <- function(colnames){
  tryCatch({
    data.frame(matrix(ncol=length(colnames), nrow=0, dimnames=list(NULL,colnames)))
  },
  error=function(e){
    cat(paste('\nSomething went wrong while creating a new Dataframe for:',colnames))
    print(e)
  })
}




write_df_to_csv <- function(df, file_path, rownames=FALSE){
  tryCatch({
      write.csv(df, file_path, row.names=rownames)
  },
  error=function(e){
    cat(paste('\nSomething went wrong while writing the file:',file_path))
    print(e)
  })
}




get_cleaned_values_from_df <- function(df){
  tryCatch({
    if (is.null(df) || (is.data.frame(df) && (nrow(df)*ncol(df)==0)) || sum(is.na(df)==nrow(df))){
      return (c())
    }else{
      # Apply some basic functions to the entire df and return the unique values without NA
      cleaned_df_values <- na.omit(unique(gsub("[[:space:]]","",trimws(as.vector(as.matrix(df))))))
      return (cleaned_df_values[!grepl("any", cleaned_df_values)])
    }
  },
  error=function(e){
    cat('\nSomething went wrong while retrieving the cleaned values of dataframe')
    print(e)
  })
}




get_num_asctb_celltypes <- function(asctb_master_table){
  tryCatch({
    
    # Initializing the dataframe of celltypes to return
    celltype_overall <- create_new_df('Cell.Types')
    
    # Get the subset dataframe of `CT/[0-9]` and `CT/[0-9]/ID` columns
    celltype_combinations <- as.data.frame(asctb_master_table[grepl("CT/[0-9]$",colnames(asctb_master_table)) | grepl("CT/[0-9]/ID$",colnames(asctb_master_table))])
    ALL_COLS <- colnames(celltype_combinations)
    n_levels <- length(ALL_COLS)/2
    
    for (i in seq(n_levels,1,-1)) {
      
      # Choose the last-level columns and check if they are completely empty
      celltypes_at_last_level <- as.data.frame(celltype_combinations[grepl(i, ALL_COLS)])
      
      if (!all(is.na(celltypes_at_last_level))){
        
        # Extract the IDs of last-level columns
        celltypes_at_last_level.ids <- as.data.frame(celltypes_at_last_level[grepl('ID',colnames(celltypes_at_last_level))])
        last_id_col <- colnames(celltypes_at_last_level.ids)
        unique_ids_at_last_level <- get_cleaned_values_from_df(celltypes_at_last_level.ids)
        
        # Append the IDs of last-level columns
        celltype_overall <- rbind(celltype_overall, setNames(as.data.frame(unique_ids_at_last_level), names(celltype_overall)))
        
        # Remove the entries which have IDs associated
        celltypes_at_last_level <- celltypes_at_last_level[is.na(celltypes_at_last_level[last_id_col]),]
        celltype_combinations <- celltype_combinations[is.na(celltype_combinations[last_id_col]),]
        
        # Extract the names of last-level columns, for the rows which didn't have any ID
        celltypes_at_last_level.names <- as.data.frame(celltypes_at_last_level[grepl('CT/[0-9]$',colnames(celltypes_at_last_level))])
        last_name_col <- colnames(celltypes_at_last_level.names)
        unique_names_at_last_level <- get_cleaned_values_from_df(celltypes_at_last_level.names)
        
        # Append the names of last-level columns
        celltype_overall <- rbind(celltype_overall, setNames(as.data.frame(unique_names_at_last_level), names(celltype_overall)))
        
        # Remove the entries which have names associated
        celltype_combinations <- celltype_combinations[is.na(celltype_combinations[last_name_col]),]
      }
    }
    
    return (as.vector(unlist(celltype_overall)))
    
  },
  error = function(e){
    cat(paste('\nSomething went wrong while generating the Celltype-vs-Counts stats for:',config$name,'\n'))
    print(e)
  })
}






process_azimuth_annotation_celltype_data <- function(body_organ, cell_hierarchy_cols) {
  tryCatch({
    # Pull the cell-annotation data from the official Azimuth Website backend-repository
    return(sapply(1:length(cell_hierarchy_cols),
                  function (i, levels) {
                    # append the individual csv file-name for the organ into the base-url of raw-github-csv files
                    ct_file <- paste0(AZIMUTH.ANNOTATION_FILES_BASE_URL, levels[i])
                    
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
    cat(paste('\nSomething went wrong while processing the cell-type annotation file for:',config$name))
    print(e)
  })
}




process_azimuth_reference <- function(config) {
  tryCatch({
    body_organ <- config$name
    cell_hierarchy_cols <- config$cell_type_columns # Still require this field for parsing the RDS data
    ref_url <- config$url # Azimuth reference RDS-file download URL
    
    # load input reference files locally if present, else download them from URL in config
    file_name <- paste0(AZIMUTH.REFERENCE_RDS_DIR, body_organ, ".Rds")
    if (!file.exists(file_name)) {
      GET(ref_url, write_disk(file_name))
    }
    ref <- readRDS(file_name)
    
    # get columns containing cell type data
    cell_types <- ref@meta.data[cell_hierarchy_cols];
    # rename columns according to ASCT+B table format
    final_column_names <- sapply(1:length(cell_hierarchy_cols), function (x) return(paste0("AS/",x)))
    names(cell_types) <- final_column_names
    
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
    cat(paste('\nSomething went wrong while processing the Azimuth reference for:',config$name))
    print(e)
  })
}




process_config_for_azimuth <- function(config) {
  tryCatch({
    # extract azimuth reference data
    reference_table <- process_azimuth_reference(config)
    
    # pull the files for cell-type ontology data from the backend Azimuth website Repo directly
    ct_ontology_tables <- process_azimuth_annotation_celltype_data(config$name, config$new_cell_type_files)
    
    # map ontology ID and LABELS to reference cell types
    # For brain the CT at L3: 'Oligo L2-6 OPALIN MAP6D1' goes missing because corresponding L4: 'exclude' is not present in annotation file.
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
    cat(paste('\nSomething went wrong while processing the config for:',config$name))
    print(e)
  })
}




get_asctb_master_table_content <- function(config){
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
    cat(paste('\nSomething went wrong while generating the ASCTB master table for:',config$name))
    print(e)
  })
}




write_asctb_structure <- function(body_organ, asctb_table) {
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
    write.table(header_rows, file =  paste0(ASCTB_TARGET_DIR, body_organ, ".csv"), 
                sep = ',', na = "", row.names = FALSE, col.names = FALSE)
    write.table(asctb_table, file = paste0(ASCTB_TARGET_DIR, body_organ, ".csv"),
                sep = ',', na = "", append = TRUE, row.names = FALSE, col.names = TRUE)
  },
  error = function(e){
    cat(paste('\nSomething went wrong while writing the ASCTB table for:',config$name))
    print(e)
  })
}
