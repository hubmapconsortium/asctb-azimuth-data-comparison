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
#                  3. Get HGNC-IDs for biomarkers:
#                      get_hgnc_id_for_biomarker() and get_all_biomarker_ids()-
#                         Gets the BG.IDs for Biomarkers whose names are not present in ASCT+B
#                  4. Write the final report for Azimuth vs ASCT+B reconciliation:
#                      create_combined_summaries()-
#                         Combines the Azimuth and ASCT+B reports into an excel spreadsheet. Currently, writes a static formula to the excel spreadsheet.
#							Could potentially be improved as a low priority item.
# AUTHOR:         Darshal Shetty/ Vikrant Deshpande/ Amber Ramesh



get_hgnc_id_for_biomarker <- function(biomarker_symbol, verbose=F){
  tryCatch({
    BASE_API_URL <- 'https://rest.genenames.org/search/symbol/'
    response <- httr::GET(httr::add_headers(`Accept` = 'application/json'),  url=paste0(BASE_API_URL, biomarker_symbol))
    
    # If API call failed return 'N/A'
    if (response$status_code!=200){
      if (verbose)  {cat('\nAPI call failed...', biomarker_symbol)}
      return ('N/A')
    }
    
    # If No match found, return 'N/A'
    response_content <- content(response)$response
    if (response_content$numFound==0){
      if (verbose)  {cat('\nNo match found', biomarker_symbol)}
      return ('N/A')
    }
    
    id <- response_content$docs[[1]]$hgnc_id
    if (verbose)  {cat('\nID found ',biomarker_symbol,' : ',id)}
    return (id)
  },
  error=function(e){
    cat('\nSomething went wrong while fetching the HGNC-ID for:',biomarker_symbol)
    print(e)
  })
}





get_all_biomarker_ids <- function(biomarkers){
  tryCatch({
    bg_ids <- c()
    for (biomarker_symbol in biomarkers){
      
      # If Biomarker name is not available in the cached dataframe, then retrieve the topmost result got from the HGNC-API and add it to cache.
      if (biomarker_symbol %in% BIOMARKER_NAME_VS_ID_MAPPING$Biomarker.Name){
        id <- BIOMARKER_NAME_VS_ID_MAPPING[BIOMARKER_NAME_VS_ID_MAPPING$Biomarker.Name==biomarker_symbol , 'Biomarker.ID']
      }else{
        id <- get_hgnc_id_for_biomarker(biomarker_symbol, verbose = T)
        if (id!='N/A'){
          BIOMARKER_NAME_VS_ID_MAPPING <<- rbind(BIOMARKER_NAME_VS_ID_MAPPING, c(biomarker_symbol, id))
        }
      }
      bg_ids <- c(bg_ids, id)
    }
    return (bg_ids)
  },
  error=function(e){
    cat('\nSomething went wrong while fetching the HGNC-ID for:',biomarker_symbol)
    print(e)
  })
}








create_new_df <- function(colnames){
  tryCatch({
    data.frame(matrix(ncol=length(colnames), nrow=0, dimnames=list(NULL,colnames)))
  },
  error=function(e){
    cat('\nSomething went wrong while creating a new Dataframe for:',colnames)
    print(e)
  })
}




write_df_to_csv <- function(df, file_path, rownames=F){
  tryCatch({
      write.csv(df, file_path, row.names=rownames, fileEncoding="UTF-8")
  },
  error=function(e){
    cat(paste('\nSomething went wrong while writing the file:',file_path))
    print(e)
  })
}




get_cleaned_values_from_df <- function(df, for_counts=T){
  tryCatch({
    if (is.null(df) || (is.data.frame(df) && (nrow(df)*ncol(df)==0)) || sum(is.na(df)==nrow(df))){
      return (c())
    }
    
    if (for_counts){
      # Removes spaces too
      cleaned_df_values <- na.omit(unique(gsub("[[:space:]]","",trimws(as.vector(as.matrix(df))))))
    }else{
      # Doesn't remove spaces
      cleaned_df_values <- na.omit(unique(trimws(as.vector(as.matrix(df)))))
    }
    
    # Remove "" blank entries, and entries which contain "any"
    cleaned_df_values <- cleaned_df_values[cleaned_df_values!=""]
    cleaned_df_values <- cleaned_df_values[!grepl("any", cleaned_df_values)]
    
    return (cleaned_df_values)
  },
  error=function(e){
    cat('\nSomething went wrong while retrieving the cleaned values of dataframe')
    print(e)
  })
}




get_hashtable_key_values <- function(df, cols, hmap.colnames, verbose=F){
  tryCatch({
    res_map <- create_new_df(hmap.colnames)
    # For each AS/1 and AS/1/ID get the unique list of ID vs Name combinations
    for (i in seq(1,length(cols),2)){
      if (verbose)  { print(paste(cols[i], cols[i+1])) }
      hmap <- aggregate(formula=as.formula(paste0('`',cols[i],'`~`',cols[i+1],'`')), data=df, FUN=unique)
      colnames(hmap) <- hmap.colnames
      res_map <- rbind(res_map, hmap)
    }
    res_map <- unique(res_map)
    return (res_map)
  },
  error=function(e){
    cat(paste('\nSomething went wrong while creating the Hashmap for ', hmap.colnames))
    print(e)
  })
}





get_cts_not_in_asctb <- function(az.hmap, cts_missing, verbose=F){
  tryCatch({
    if (verbose)  {print(paste('Unlisting the Cell-names since some Azimuth CTs have multiple names for same ID.'))}
    az.hmap$Cell.Names <- lapply(az.hmap$Cell.Names, `[[`, 1)
    az.cts_not_in_asctb <- left_join(cts_missing, az.hmap)
    
    if (verbose)  {print(paste('Replacing NA/NULL with "N/A" after left-joining with Azimuth'))}
    az.cts_not_in_asctb$Cell.Names <- unlist(replace_na(az.cts_not_in_asctb$Cell.Names, 'N/A'))
    
    if (verbose)  {print(paste('Remove CTs with junk characters not like %CT:% and then perform final aggregation to choose one name for one CT-ID.'))}
    az.cts_not_in_asctb <- az.cts_not_in_asctb[grepl('CL:',az.cts_not_in_asctb$Cell.IDs) ,]
    az.cts_not_in_asctb <- aggregate(formula=as.formula('Cell.Names~Cell.IDs'), data=az.cts_not_in_asctb, FUN=max)
    return (az.cts_not_in_asctb)
  },
  error=function(e){
    cat(paste('\nSomething went wrong with left-join for missing ASCTB-CTs with Azimuth Hashmap of CTs'))
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
                  simplify = F))
  },
  error = function(e){
    cat('\nSomething went wrong while processing the cell-type annotation file for:',config$name)
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
    cat('\nSomething went wrong while processing the Azimuth reference for:',config$name)
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
    merged_data <- Reduce(left_join, ct_ontology_tables, reference_table)
    
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
    cat('\nSomething went wrong while processing the config for:',config$name)
    print(e)
  })
}






get_asctb_master_table_content <- function(config){
  # Read the google-sheet tab url and return a dataframe after replacing the metainformation in top 10 rows. If URL not available return NA.
  tryCatch({
    
    url <- config$asctb_master_url
    if (url=="NA")          { return (NA) }
    
    asctb.master.data <- gsheet2tbl(url)
    
    # Remove out the top 10 meta info rows
    asctb.master.data <- asctb.master.data[10:nrow(asctb.master.data),]
    
    # Remove out the top row which was just the colnames
    colnames(asctb.master.data) <- asctb.master.data[1,]
    asctb.master.data <- as.data.frame(asctb.master.data[2:nrow(asctb.master.data),])
    
    file_path <- paste0(STAGING_DIR,config$name,'_asctb_master.csv')
    write_df_to_csv(asctb.master.data, file_path)
    
    return (file_path)
    
  },
  error=function(e){
    cat('\nSomething went wrong while generating the ASCTB master table for:',config$name)
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
                          byrow = T)
    write.table(header_rows, file =  paste0(ASCTB_TARGET_DIR, body_organ, ".csv"), 
                sep = ',', na = "", row.names = F, col.names = F)
    write.table(asctb_table, file = paste0(ASCTB_TARGET_DIR, body_organ, ".csv"),
                sep = ',', na = "", append = T, row.names = F, col.names = T)
  },
  error = function(e){
    cat('\nSomething went wrong while writing the ASCTB table for:',config$name)
    print(e)
  })
}






create_combined_summaries <- function(asctb_organ_stats, azimuth_organ_stats, verbose=F){
  tryCatch({
    
    azimuth_organ_stats <- azimuth_organ_stats[order(azimuth_organ_stats$Organ),]
    asctb_organ_stats <- asctb_organ_stats[order(asctb_organ_stats$Organ),]
    
    combined_report.cols_ordered <- c("Organ", "AZ.Annotation.Levels", "AZ.Unique.CTs", "AZ.Unique.CT.IDs", "ASCTB.Unique.CTs", "ASCTB.Unique.CT.IDs", "Matching.CT.IDs", "CTwID.Missing.in.ASCTB", 
                                      "AZ.Total.Cells", "AZ.Unique.BGs", "ASCTB.Unique.BGs", "Matching.BGs", "BGwID.Missing.in.ASCTB",  "Raw.Organ.Name")
    combined.azimuth_vs_asctb <- left_join(asctb_organ_stats, azimuth_organ_stats, by="Organ")
    combined.azimuth_vs_asctb <- combined.azimuth_vs_asctb[,combined_report.cols_ordered]
    
    
    # Add a final row for Totals
    res_row <- c()
    for (col in colnames(combined.azimuth_vs_asctb)){
      if (grepl('Organ',col))       { res_row <- c(res_row, 'Totals:') }
      else if (grepl('Link', col))  { res_row <- c(res_row, '') }
      else                          { res_row <- c( res_row, sum(as.numeric(unlist(combined.azimuth_vs_asctb[col]))) ) }
    }
    combined.azimuth_vs_asctb <- rbind(combined.azimuth_vs_asctb, res_row)
    
    # Attempting to create hyperlinks in excel file. Trim the length of each filename since sheet-name can have <=31 characters
    combined.azimuth_vs_asctb['Shortcut for CTs not in ASCT+B'] <- c(as.vector(sapply( as.vector(unlist(combined.azimuth_vs_asctb['Raw.Organ.Name'])),
                                                                  function(organ) {
                                                                    if (organ=='Totals:' || combined.azimuth_vs_asctb[combined.azimuth_vs_asctb$Raw.Organ.Name==organ,]$CTwID.Missing.in.ASCTB=="0"){
                                                                      return ('')
                                                                      }else{
                                                                        sheet_name <- substr(paste0(organ,'.cts_not_in_asctb'), 1, 27)
                                                                        return (paste0('=HYPERLINK("#', sheet_name, '!A1", "Missing CTs")'))
                                                                        }
                                                                    })))
     
    # combined.azimuth_vs_asctb['Shortcut for BGs not in ASCT+B'] <- c(as.vector(sapply( as.vector(unlist(combined.azimuth_vs_asctb['Raw.Organ.Name'])),
    #                                                               function(organ) {
    #                                                                 if (organ=='Totals:' || combined.azimuth_vs_asctb[combined.azimuth_vs_asctb$Raw.Organ.Name==organ,]$BGwID.Missing.in.ASCTB=="0"){
    #                                                                   return ('')
    #                                                                 }else{
    #                                                                   sheet_name <- substr(paste0(organ,'.bgs_not_in_asctb'), 1, 27)
    #                                                                   return (paste0('=HYPERLINK("#', sheet_name, '!A1", "Missing BGs")'))
    #                                                                 }
    #                                                               })))
     
    
    
    
    # Read all files in the staging directory that contain Azimuth minus ASCTB information
    files <- list.files(STAGING_DIR)
    files <- sort(files[grepl("cts_not_in_asctb.csv|bgs_not_in_asctb.csv", files)])
    lst <- list()
    lst[['Azimuth_vs_ASCTB']] <- combined.azimuth_vs_asctb[, names(combined.azimuth_vs_asctb)!=c("Raw.Organ.Name")]
    
    # Trim the length of each filename since sheet-name can have <=31 characters
    for (i in 1:length(files)){
      if (verbose)  {print(files[i])}
      df <- read.csv(paste0(STAGING_DIR,files[i]))
      sheet_name <- sub( '.csv', '', files[i])
      sheet_name <- substr(sheet_name, 1, 27)
      
      sorted_cols <- sort(colnames(df))
      df <- as.data.frame(df[, sorted_cols])
      colnames(df) <- sorted_cols
      lst[[sheet_name]] <- df
    }
    
    
    write_df_to_csv(azimuth_organ_stats, paste0(SUMMARIES_DIR,'Azimuth.All_organs.stats.csv'))
    write_df_to_csv(asctb_organ_stats, paste0(SUMMARIES_DIR,'ASCTB.All_organs.stats.csv'))
    write.xlsx(lst, file=paste0(SUMMARIES_DIR, 'Azimuth_vs_ASCTB.summaries.xlsx'), overwrite=T)
    
    
    # Also write the latest cache to the Biomarker-name-vs-id csv
    BIOMARKER_NAME_VS_ID_MAPPING <<- unique(BIOMARKER_NAME_VS_ID_MAPPING)
    write_df_to_csv(BIOMARKER_NAME_VS_ID_MAPPING, BIOMARKER_NAME_VS_ID_CACHE)
    
},
error = function(e){
  cat('\nSomething went wrong while writing the combined summary file for Azimuth vs ASCT+B.')
  print(e)
})
}
