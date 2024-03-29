# DOCSTRING:      This script has summary computation functions for:
#                   1. The Azimuth Reference data that was converted into ASCTB compatible structure.
#                   2. The ASCT+B Master datasets.
#                   3. All Cell-Types within a specific Azimuth reference.
#                   4. Getting celltype counts using an iterative algorithm that traverses the dataset from right to left.
#						This was necessary to match counts against the CCF-reporter counts which parses JSON, but here we have CSV formatted tabular data.
# AUTHOR:         Vikrant Deshpande / Amber Ramesh





generate_azimuth_celltype_cnt_summary <- function(asct_table, organ){
  tryCatch({
    
    # Structure for stats table of organ level CellType vs count
    celltype_vs_counts_cols <- c('Cell.Type', 'Cell.Type.ID', 'Num.Cells', 'Annotation.Level')
    celltype_vs_counts <- create_new_df(celltype_vs_counts_cols)
    
    # Extract the available CellType columns; Azimuth dataset has CT-columns named as 'AS/'...
    celltype_combinations <- asct_table[grepl("AS/[0-9]$|AS/[0-9]/ID$|COUNT$",colnames(asct_table))]
    ALL_COLS <- colnames(celltype_combinations)
    COUNT_COL <- ALL_COLS[grepl("COUNT",ALL_COLS)]
    ID_COLS <- ALL_COLS[grepl("ID", ALL_COLS)]
    ALL_COLS <- ALL_COLS[ALL_COLS!=COUNT_COL]
    
    
    # Edge case: Spleen doesn't have any ID columns, Fetal_Development doesn't have any Cell-Type IDs
    if (length(ID_COLS)==0){
      ID_COLS <- sapply(1:length(ALL_COLS), function (x) return(paste0("AS/",x,"/ID")))
      celltype_combinations[ID_COLS] <- "NA"
      ALL_COLS <- c(ALL_COLS, ID_COLS)
    }else if(all(is.na(celltype_combinations[, ID_COLS]))){
      celltype_combinations <- celltype_combinations[, !grepl("ID", ALL_COLS)]
      celltype_combinations[ID_COLS] <- "NA"
    }
    ALL_COLS <- sort(ALL_COLS)
    celltype_combinations <- celltype_combinations[c(ALL_COLS, COUNT_COL)]
    
    for (i in (length(ALL_COLS)/2):1) {
      
      # Choose the (keyCols, COUNT) combo -> aggregate it -> merge it with the CellType stats table created above
      celltypes_at_last_level <- as.data.frame(celltype_combinations[ , c(ALL_COLS, COUNT_COL)])
      
      # Group by all the key-columns and return the value in the COUNT_COL
      group_by_formula <- as.formula(paste0( paste0('`',COUNT_COL,'`~'), paste0(sapply(ALL_COLS, function(col) paste0('`',col,'`')), collapse="+")))
      agg_celltypes_at_last_level <- aggregate( group_by_formula, data=celltypes_at_last_level, FUN = sum )
      
      ncols <- ncol(agg_celltypes_at_last_level)
      celltype_counts_last_level <- agg_celltypes_at_last_level[ , c(ncols-2, ncols-1, ncols)]
      
      agg_celltypes <- celltype_counts_last_level[!is.na(celltype_counts_last_level[,1]) & !is.null(celltype_counts_last_level[,1]), ]
      agg_celltypes['Annotation.Level'] <- i
      agg_celltypes <- agg_celltypes[!trimws(agg_celltypes[,1]) %in% trimws(celltype_vs_counts[,1]),]
      
      # Append the output of current key-columns aggregation into the stats table
      celltype_vs_counts <- rbind(celltype_vs_counts, setNames(agg_celltypes, names(celltype_vs_counts)))
      
      # Remove the last column and move on to aggregating the next subset of columns
      ALL_COLS <- ALL_COLS[!grepl(i,ALL_COLS)]
    }
    
    # Perform a final aggregation since Brain's annotation-files have L4 level issues.
    celltype_vs_counts <- aggregate( as.formula(paste0('`Num.Cells`~`Cell.Type`+`Cell.Type.ID`+`Annotation.Level`')), data=celltype_vs_counts, FUN = sum )
    
    # Reformat the stats dataframe and write it to a csv
    celltype_vs_counts['Organ'] <- organ
    celltype_vs_counts <- celltype_vs_counts[order(celltype_vs_counts$Cell.Type), c('Organ', 'Cell.Type', 'Cell.Type.ID', 'Annotation.Level', 'Num.Cells')]
    write_df_to_csv(celltype_vs_counts, paste0(SUMMARIES_DIR, organ, '.celltype_stats.csv'))
  },
  error = function(e){
    cat('\nSomething went wrong while generating the Celltype-vs-Counts stats for:',config$name,'\n')
    print(e)
  })
}




process_azimuth_ref_dataset_summary <- function(config, asct_table, azimuth_organ_stats){
  tryCatch({
    # Wrangle the Azimuth dataset to derive summary stats
    asctb.formatted_azimuth_columns <- colnames(asct_table)
    
    # C1: Organ Name
    organ.1 <- config$asctb_name
    
    # Generate a separate summary table: celltype, count(*) for this organ
    generate_azimuth_celltype_cnt_summary(asct_table, config$name)
    
    # C2: Get the union of all "CT" columns in the ASCTB organ.csv file
    entire_set_of_cell_types <- asct_table[grepl("AS/[0-9]$",asctb.formatted_azimuth_columns)]
    cleaned_set_of_cell_types <- get_cleaned_values_from_df(entire_set_of_cell_types)
    n_unique_cell_types.2 <- length(cleaned_set_of_cell_types)
    
    # C3: Get the union of all "ID" columns in the ASCTB organ.csv file
    entire_set_of_ct_ontology_ids <- asct_table[grepl("AS",asctb.formatted_azimuth_columns) & grepl("ID",asctb.formatted_azimuth_columns)]
    cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(entire_set_of_ct_ontology_ids)
    n_unique_ct_ontology_ids.3 <- length(cleaned_set_of_ct_ontology_ids)
    
    # C4: Get the sum of all counts that were already generated using the groupby formula used in process_config_for_azimuth()
    n_total_cell_types.4 <- sum(asct_table[grep('COUNT',asctb.formatted_azimuth_columns)])
    
    # C5: Get the number of annotations for an organ directly from the json config file
    n_annotation_levels.5 <- length(config$new_cell_type_files)
    
    # C6: Get the global variable containing all biomarkers for an organ, that we set during file-processing
    azimuth.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(azimuth.entire_set_of_biomarkers)
    n_unique_biomarkers.6 <- length(azimuth.cleaned_set_of_biomarkers)
    
    # Append (C1, C2, ... C6) to the existing azimuth_organ_stats global var
    azimuth_organ_stats <<- rbind(azimuth_organ_stats, c( organ.1, n_unique_cell_types.2, n_unique_ct_ontology_ids.3, 
                                                          n_total_cell_types.4, n_annotation_levels.5, n_unique_biomarkers.6, config$name))
    colnames(azimuth_organ_stats) <<- azimuth_organ_stats_cols
  },
  error = function(e){
    cat('\nSomething went wrong while generating the Azimuth organ-level stats for:',config$name,"\n")
    print(e)
  })
}






get_num_asctb_celltypes <- function(asctb_master_table, organ, verbose=F){
  tryCatch({
    
    # Get the subset dataframe of `CT/[0-9]` and `CT/[0-9]/ID` columns.
    celltype_combinations <- as.data.frame(asctb_master_table[grepl("CT/[0-9]$|CT/[0-9]/ID$",colnames(asctb_master_table))])
    ALL_COLS <- colnames(celltype_combinations)
    
    if (organ %in% c("Brain")){
      return (get_cleaned_values_from_df(asctb_master_table[ALL_COLS]))
    }
    
    # Initializing the dataframe of celltypes to return
    hmap.colnames <- c("Cell.IDs", "Cell.Names")
    celltype_overall <- create_new_df(hmap.colnames)
    
    for (i in seq(1,length(ALL_COLS),2)) {
      if (verbose)  { print(paste(ALL_COLS[i], ALL_COLS[i+1])) }
      
      hmap <- get_hashtable_key_values(asctb_master_table, c(ALL_COLS[i], ALL_COLS[i+1]), hmap.colnames)
      celltype_overall <- rbind(celltype_overall, hmap)
      cts.without.ids <- as.data.frame( get_cleaned_values_from_df(asctb_master_table[is.na(asctb_master_table[ALL_COLS[i+1]]), ALL_COLS[i]]) )
      if (nrow(cts.without.ids)>0 && !all(is.na(cts.without.ids)) && !all(is.null(cts.without.ids))){
        if (verbose)  { print('Appending the CTs which have no names...') }
        cts.without.ids <- data.frame(rep("",nrow(cts.without.ids)), cts.without.ids)
        colnames(cts.without.ids) <- hmap.colnames
        celltype_overall <- rbind(celltype_overall, cts.without.ids)
      }
    }
    
    celltype_overall <- unique(celltype_overall)
    return (get_cleaned_values_from_df(unlist(celltype_overall$Cell.Names)))
  },
  error = function(e){
    if (verbose)  {traceback()}
    cat('\nSomething went wrong while getting the Celltypes for:',config$name,'\n')
    print(e)
  })
}





process_asctb_master_dataset_summary <- function(config, file_path, asct_table_derived_from_azimuth){
  tryCatch({
    
    # Wrangle the ASCT+B dataset to derive summary stats, or just add a dummy entry when no ASCTB-Master table
    if (is.na(file_path) || !file.exists(file_path)){
      asctb_organ_stats <<- rbind(asctb_organ_stats, c(config$asctb_name, rep(0,length(asctb_organ_stats_cols)-1)))
      colnames(asctb_organ_stats) <<- asctb_organ_stats_cols
      return(paste0("\nAppended default values for ",config$asctb_name,"."))
    }
    


    asctb_master_table <- as.data.frame(read.csv(file_path, na.string=c("NA", "NULL"), encoding="UTF-8"))
    colnames(asctb_master_table) <- gsub('\\.', '/', colnames(asctb_master_table))
    asctb.master_columns <- colnames(asctb_master_table)
    
    # C1: Organ Name
    organ.1 <- config$asctb_name
    
    # C2: Get the union of all "CT" columns in the ASCTB organ.csv file
    cat('\nGenerating the set of cell-types now')
    asctb.cell_types <- get_num_asctb_celltypes(asctb_master_table, organ.1)
    cat('\nNow, cell-types for ', organ.1,' are: ')
    print(asctb.cell_types)
    n_unique_cell_types.2 <- length(asctb.cell_types)
    
    
    # C3: Get the union of all "ID" columns in the ASCTB organ.csv file. Kidney has extra 'CT' column that messes up the counts
    asctb.entire_set_of_ct_ontology_ids <- asctb_master_table[(grepl("CT",asctb.master_columns)) & grepl("ID",asctb.master_columns) & !grepl("AS",asctb.master_columns)]
    asctb.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(asctb.entire_set_of_ct_ontology_ids)
    n_unique_ct_ontology_ids.3 <- length(asctb.cleaned_set_of_ct_ontology_ids)
    
    
    # C4: Get the intersection of CT_Ontology_IDs in the ASCTB organ.csv file vs Azimuth organ.csv file.
    azimuth_colnames <- colnames(asct_table_derived_from_azimuth)
    
    
    # Edge Case: ASCTB Brain/Spleen doesn't have CT_Ontology_IDs. Use Cell-Type Names for finding # matching between ASCTB and Azimuth.
    if (organ.1 %in% c('Spleen','Brain')){
      asctb.cleaned_set_of_cell_types <- get_cleaned_values_from_df(asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns)])
      
      azimuth.entire_set_of_cell_types <- asct_table_derived_from_azimuth[grepl("AS/[0-9]$",azimuth_colnames)]
      azimuth.cleaned_set_of_cell_types <- get_cleaned_values_from_df(azimuth.entire_set_of_cell_types)
      
      if (organ.1=='Brain'){
        
        n_matching_ct_ontology_ids.4 <- length(intersect(asctb.cleaned_set_of_cell_types, azimuth.cleaned_set_of_cell_types))
        azimuth_cts_missing_in_asctb <- as.data.frame(setdiff( get_cleaned_values_from_df(azimuth.entire_set_of_cell_types,F) ,
                                                               get_cleaned_values_from_df(asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns)],F)))
        colnames(azimuth_cts_missing_in_asctb) <- c("Cell.Names")
        cols <- sort(azimuth_colnames[grepl("AS/[0-9]/ID$|AS/[0-9]$",azimuth_colnames)])
        hmap.cols <- c("Cell.IDs", "Cell.Names")
        hmap <- get_hashtable_key_values(asct_table_derived_from_azimuth, cols, hmap.cols)
        hmap <- as.data.frame(hmap %>% unnest(Cell.Names, keep_empty=T))
        azimuth_cts_not_in_asctb <- unique(left_join(azimuth_cts_missing_in_asctb, hmap))
      }else{
        n_matching_ct_ontology_ids.4 <- length(intersect( sapply(get_cleaned_values_from_df(azimuth.entire_set_of_cell_types,F), tolower) , 
                                                                  sapply(get_cleaned_values_from_df(asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns)],F), tolower)))
        azimuth_cts_not_in_asctb <- as.data.frame(setdiff( sapply(get_cleaned_values_from_df(azimuth.entire_set_of_cell_types,F), tolower) , 
                                                               sapply(get_cleaned_values_from_df(asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns)],F), tolower)))
        colnames(azimuth_cts_not_in_asctb) <- c("Cell.Names")
      }
      
    }else{
      azimuth.entire_set_of_ct_ontology_ids <- asct_table_derived_from_azimuth[grepl("AS/[0-9]/ID$",azimuth_colnames)]
      azimuth.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(azimuth.entire_set_of_ct_ontology_ids)
      n_matching_ct_ontology_ids.4 <- length(intersect(asctb.cleaned_set_of_ct_ontology_ids, azimuth.cleaned_set_of_ct_ontology_ids))
      
      azimuth.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(azimuth.entire_set_of_ct_ontology_ids, for_counts=F)
      azimuth_cts_missing_in_asctb <- as.data.frame(setdiff(azimuth.cleaned_set_of_ct_ontology_ids , asctb.cleaned_set_of_ct_ontology_ids))
      colnames(azimuth_cts_missing_in_asctb) <- c("Cell.IDs")
      
      cols <- sort(azimuth_colnames[grepl("AS/[0-9]/ID$|AS/[0-9]$",azimuth_colnames)])
      hmap.cols <- c("Cell.IDs", "Cell.Names")
      hmap <- get_hashtable_key_values(asct_table_derived_from_azimuth, cols, hmap.cols)
      azimuth_cts_not_in_asctb <- get_cts_not_in_asctb(az.hmap=hmap, cts_missing=azimuth_cts_missing_in_asctb)
    }
    
    n_missing_ct_ontology_ids.4.2 <- nrow(azimuth_cts_not_in_asctb)
    write_df_to_csv(df=azimuth_cts_not_in_asctb, file_path=paste0(STAGING_DIR, config$name, '.cts_not_in_asctb.csv'))
    
    
    
    # C5: Get the union of all "BG" columns in the ASCTB organ.csv file
    asctb.entire_set_of_biomarkers <<- asctb_master_table[(grepl("BG",asctb.master_columns)) & !(grepl("ID|COUNT|LABEL",asctb.master_columns))]
    asctb.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(asctb.entire_set_of_biomarkers)
    asctb.cleaned_biomarker_ids <- get_cleaned_values_from_df(asctb_master_table[(grepl("BG",asctb.master_columns) & grepl("ID",asctb.master_columns)) & !(grepl("COUNT|LABEL",asctb.master_columns))])
    n_unique_biomarkers.5 <- length(asctb.cleaned_set_of_biomarkers)
    
    asctb.entire_set_of_biomarkers_prots <- asctb_master_table[(grepl("BP",asctb.master_columns)) & !(grepl("ID|COUNT|LABEL",asctb.master_columns))]
    asctb.cleaned_set_of_biomarkers_prots <- get_cleaned_values_from_df(asctb.entire_set_of_biomarkers_prots)
    n_unique_biomarkers.5.2 <- length(asctb.cleaned_set_of_biomarkers_prots)
    
    
    # C6: Get the intersection of Biomarkers in the ASCTB organ.csv file vs Azimuth organ.csv file
    azimuth.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(azimuth.entire_set_of_biomarkers)
    n_matching_biomarkers.6 <- length(intersect(azimuth.cleaned_set_of_biomarkers, asctb.cleaned_set_of_biomarkers))
    azimuth_bgs_not_in_asctb <- setdiff(azimuth.cleaned_set_of_biomarkers, c(asctb.cleaned_set_of_biomarkers, asctb.cleaned_set_of_biomarkers_prots))
    
    n_missing_biomarkers.6.2 <- 0
    
    if (!is.null(azimuth_bgs_not_in_asctb) && !all(is.na(azimuth_bgs_not_in_asctb))){
      print(paste0('Retrieving the HGNC-IDs for ',length(azimuth_bgs_not_in_asctb),' Azimuth Biomarker-names not present in ASCTB...'))
      bg_ids <- get_all_biomarker_ids(azimuth_bgs_not_in_asctb)
      azimuth_bgs_not_in_asctb <- data.frame(bg_ids, azimuth_bgs_not_in_asctb)
      colnames(azimuth_bgs_not_in_asctb) <- c('Biomarker.IDs','Biomarker.Names')
      
      # After comparing Azimuth-BG names with the ASCTB-BG names, now remove the Azimuth-BG IDs present in ASCTB
      azimuth_bg_ids_in_asctb <- azimuth_bgs_not_in_asctb$Biomarker.IDs %in% asctb.cleaned_biomarker_ids
      n_matching_biomarkers.6 <- n_matching_biomarkers.6 + nrow(azimuth_bgs_not_in_asctb[azimuth_bg_ids_in_asctb ,])
      azimuth_bgs_not_in_asctb <- azimuth_bgs_not_in_asctb[!azimuth_bg_ids_in_asctb ,]
      n_missing_biomarkers.6.2 <- nrow(azimuth_bgs_not_in_asctb)
      
      write_df_to_csv(df=azimuth_bgs_not_in_asctb, file_path=paste0(STAGING_DIR, config$name, '.bgs_not_in_asctb.csv'))
    }
    
    # Append (C1, C2, ... C6) to the existing asctb_organ_stats global var for intersection-stats.
    asctb_organ_stats <<- rbind(asctb_organ_stats, c( organ.1, n_unique_cell_types.2, n_unique_ct_ontology_ids.3, n_matching_ct_ontology_ids.4, n_missing_ct_ontology_ids.4.2, 
                                                      n_unique_biomarkers.5, n_unique_biomarkers.5.2, n_matching_biomarkers.6, n_missing_biomarkers.6.2 ))
    colnames(asctb_organ_stats) <<- asctb_organ_stats_cols
    
    return(paste0("Computed statistics for ",config$asctb_name,"."))
  },
  error = function(e){
    cat('\nSomething went wrong while generating the ASCTB organ-level stats for:',config$name,'\n')
    print(e)
  })
}
