# DOCSTRING:      This script has summary computation functions for both: 
#                   1. The Azimuth Reference data that was converted into ASCTB compatible structure.
#                   2. The ASCT+B Master datasets
# AUTHOR:         Vikrant Deshpande / Amber Ramesh



generate_azimuth_celltype_cnt_summary <- function(asct_table, body_organ){
  tryCatch({
    
    # Structure for stats table of organ level CellType vs count
    celltype_vs_counts_cols <- c('Cell.Types', 'Num.Cells')
    celltype_vs_counts <- create_new_df(celltype_vs_counts_cols)
    
    # Convert the available CellType columns
    celltype_combinations <- asct_table[grepl("AS/[0-9]$",colnames(asct_table)) | grepl("COUNT$",colnames(asct_table))]
    ALL_COLS <- colnames(celltype_combinations)
    COUNT_COL <- ALL_COLS[grepl("COUNT",ALL_COLS)]
    ALL_COLS <- ALL_COLS[ALL_COLS!=COUNT_COL]
    
    for (column in rev(ALL_COLS)) {
      # Choose the (keyCols, COUNT) combo -> aggregate it -> merge it with the CellType stats table created above
      celltypes_at_last_level <- celltype_combinations[ , c(ALL_COLS, COUNT_COL)]
      
      # Group by all the key-columns and return the value in the COUNT_COL
      group_by_formula <- as.formula(paste0( paste0('`',COUNT_COL,'`~'), paste0(sapply(ALL_COLS, function(col) paste0('`',col,'`')), collapse="+")))
      agg_celltypes_at_last_level <- aggregate( group_by_formula, data=celltypes_at_last_level, FUN = sum )
      
      n <- ncol(agg_celltypes_at_last_level)
      celltype_counts_last_level <- agg_celltypes_at_last_level[ , c(n-1, n)]
      agg_celltypes <- celltype_counts_last_level[!is.na(celltype_counts_last_level[,1]) & !is.null(celltype_counts_last_level[,1]), ]
      agg_celltypes <- agg_celltypes[!trimws(agg_celltypes[,1]) %in% trimws(celltype_vs_counts[,1]),]
      # Append the output of current key-columns aggregation into the stats table
      celltype_vs_counts <- rbind(celltype_vs_counts, setNames(agg_celltypes, names(celltype_vs_counts)))
      
      # Remove the last column and move on to aggregating the next subset of columns
      ALL_COLS <- ALL_COLS[ALL_COLS!=column]
    }
    
    # Perform a final aggregation since Brain's annotation-files have L4 level issues.
    celltype_vs_counts <- aggregate( as.formula(paste0('`Num.Cells`~`Cell.Types`')), data=celltype_vs_counts, FUN = sum )
    
    # Reformat the stats dataframe and write it to a csv
    celltype_vs_counts <- setNames(celltype_vs_counts, celltype_vs_counts_cols)
    celltype_vs_counts['Organ'] <- body_organ
    celltype_vs_counts <- celltype_vs_counts[, c("Organ",celltype_vs_counts_cols)]
    write_df_to_csv(celltype_vs_counts, paste0(SUMMARIES_DIR, body_organ, ".celltype_stats.csv"))
  },
  error = function(e){
    cat(paste('\nSomething went wrong while generating the Celltype-vs-Counts stats for:',config$name,'\n'))
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
                                                          n_total_cell_types.4, n_annotation_levels.5, n_unique_biomarkers.6))
    colnames(azimuth_organ_stats) <<- azimuth_organ_stats_cols
  },
  error = function(e){
    cat(paste('\nSomething went wrong while generating the Azimuth organ-level stats for:',config$name,"\n"))
    print(e)
  })
}




process_asctb_master_dataset_summary <- function(config, asctb_master_table, asct_table_derived_from_azimuth, compute_intersection_stats=TRUE){
  tryCatch({
    # Wrangle the ASCT+B dataset to derive summary stats, or just add a dummy entry when no ASCTB-Master table
    if(is.na.data.frame(asctb_master_table)){
      asctb_organ_stats <<- rbind(asctb_organ_stats, c(config$asctb_name, rep(0,length(asctb_organ_stats_cols)-1)))
      colnames(asctb_organ_stats) <<- asctb_organ_stats_cols
      return(paste0("\nAppended default values for ",config$asctb_name,"."))
    }
    
    asctb.master_columns <- colnames(asctb_master_table)
    
    # C1: Organ Name
    organ.1 <- config$asctb_name
    
    # C2: Get the union of all "CT" columns in the ASCTB organ.csv file
    asctb.entire_set_anatomical_structs <- asctb_master_table[grepl("AS/[0-9]$",asctb.master_columns)]
    asctb.cleaned_set_of_anatomical_structs <- get_cleaned_values_from_df(asctb.entire_set_anatomical_structs)
    n_unique_anatomical_structs.2 <- length(asctb.cleaned_set_of_anatomical_structs)
    
    # C3: Get the union of all "CT" columns in the ASCTB organ.csv file
    # Previous logic was incorrect.
    # asctb.entire_set_of_cell_types <- asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns)]
    # asctb.cleaned_set_of_cell_types <- get_cleaned_values_from_df(asctb.entire_set_of_cell_types)
    print('Generating the set of cell-types now')
    asctb.cleaned_set_of_cell_types <- get_num_asctb_celltypes(asctb_master_table)
    n_unique_cell_types.3 <- length(asctb.cleaned_set_of_cell_types)
    
    
    # C4: Get the union of all "ID" columns in the ASCTB organ.csv file. Kidney has extra 'CT' column that messes up the counts
    asctb.entire_set_of_ct_ontology_ids <- asctb_master_table[(grepl("CT",asctb.master_columns)) & grepl("ID",asctb.master_columns)]
    asctb.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(asctb.entire_set_of_ct_ontology_ids)
    n_unique_ct_ontology_ids.4 <- length(asctb.cleaned_set_of_ct_ontology_ids)
    
    
    if (compute_intersection_stats){
      # C5: Get the intersection of CT_Ontology_IDs in the ASCTB organ.csv file vs Azimuth organ.csv file.
      # Edge Case: Brain/Spleen doesn't have CT_Ontology_IDs. Use Cell-Types itself for finding # matching between ASCTB and Azimuth.
      azimuth_colnames <-  colnames(asct_table_derived_from_azimuth)
      if (organ.1 %in% c('Spleen','Brain')){
        azimuth.entire_set_of_cell_types <- asct_table_derived_from_azimuth[grepl("AS/[0-9]$",azimuth_colnames)]
        azimuth.cleaned_set_of_cell_types <- get_cleaned_values_from_df(azimuth.entire_set_of_cell_types)
        n_matching_ct_ontology_ids.5 <- length(intersect(asctb.cleaned_set_of_cell_types, azimuth.cleaned_set_of_cell_types))
      }else{
        azimuth.entire_set_of_ct_ontology_ids <- asct_table_derived_from_azimuth[grepl("AS",azimuth_colnames) & grepl("ID",azimuth_colnames)]
        azimuth.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(azimuth.entire_set_of_ct_ontology_ids)
        n_matching_ct_ontology_ids.5 <- length(intersect(asctb.cleaned_set_of_ct_ontology_ids, azimuth.cleaned_set_of_ct_ontology_ids))
      }
    }
    
    
    # C6: Get the union of all "BG" columns in the ASCTB organ.csv file
    asctb.entire_set_of_biomarkers <<- asctb_master_table[(grepl("BG",asctb.master_columns)) & !(grepl("ID",asctb.master_columns) | grepl("COUNT",asctb.master_columns) | grepl("LABEL",asctb.master_columns))]
    asctb.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(asctb.entire_set_of_biomarkers)
    n_unique_biomarkers.6 <- length(asctb.cleaned_set_of_biomarkers)
    
    
    if (compute_intersection_stats){
      # C7: Get the intersection of Biomarkers in the ASCTB organ.csv file vs Azimuth organ.csv file
      azimuth.cleaned_set_of_biomarkers <- get_cleaned_values_from_df(azimuth.entire_set_of_biomarkers)
      n_matching_biomarkers.7 <- length(intersect(azimuth.cleaned_set_of_biomarkers, asctb.cleaned_set_of_biomarkers))
    }
    
    
    # Append (C1, C2, ... C6) to the existing asctb_organ_stats global var if intersection-stats to be included. Otherwise append only unique counts for summary.
    if (compute_intersection_stats){
      asctb_organ_stats <<- rbind(asctb_organ_stats, c( organ.1, n_unique_anatomical_structs.2, n_unique_cell_types.3, n_unique_ct_ontology_ids.4,
                                                      n_matching_ct_ontology_ids.5, n_unique_biomarkers.6, n_matching_biomarkers.7 ))
      colnames(asctb_organ_stats) <<- asctb_organ_stats_cols
    }else{
      asctb_all_organ_summaries <<- rbind(asctb_all_organ_summaries, c( organ.1, n_unique_anatomical_structs.2, n_unique_cell_types.3, n_unique_ct_ontology_ids.4, n_unique_biomarkers.6 ))
      colnames(asctb_organ_stats) <<- asctb_all_organ_summaries_cols
    }
    return(paste0("Computed statistics for ",config$asctb_name,"."))
  },
  error = function(e){
    cat(paste('\nSomething went wrong while generating the ASCTB organ-level stats for:',config$name,'\n'))
    print(e)
  })
}
