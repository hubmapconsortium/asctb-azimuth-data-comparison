# DOCSTRING:      This script has summary computation functions for both: 
#                   1. The Azimuth Reference data that was converted into ASCTB compatible structure.
#                   2. The ASCT+B Master datasets
# AUTHOR:         Vikrant Deshpande / Amber Ramesh




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
    cat(paste('\nSomething went wrong while generating the Celltype-vs-Counts summary for:',config$name,'\n'))
    print(e)
  })
}




process_azimuth_ref_dataset_summary <- function(config, asct_table, azimuth_organ_stats){
  tryCatch({
    # Wrangle the Azimuth dataset to derive summary stats
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
    cat(paste('\nSomething went wrong while generating the Azimuth organ-level summary stats for:',config$name,"\n"))
    print(e)
  })
}




process_asctb_master_dataset_summary <- function(config, asctb_master_table, asctb_organ_stats, asct_table_derived_from_azimuth){
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
    asctb.entire_set_of_cell_types <- asctb_master_table[grepl("CT/[0-9]$",asctb.master_columns)]
    asctb.cleaned_set_of_cell_types <- get_cleaned_values_from_df(asctb.entire_set_of_cell_types)
    n_unique_cell_types.2 <- length(asctb.cleaned_set_of_cell_types)
    
    
    # C3: Get the union of all "ID" columns in the ASCTB organ.csv file. Kidney has extra 'CT' column that messes up the counts
    asctb.entire_set_of_ct_ontology_ids <- asctb_master_table[(grepl("CT",asctb.master_columns)) & grepl("ID",asctb.master_columns)]
    asctb.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(asctb.entire_set_of_ct_ontology_ids)
    n_unique_ct_ontology_ids.3 <- length(asctb.cleaned_set_of_ct_ontology_ids)
    
    
    # C4: Get the intersection of CT_Ontology_IDs in the ASCTB organ.csv file vs Azimuth organ.csv file.
      # Edge Case: Brain doesn't have any CT_Ontology_IDs. Use Cell-Types itself for finding # matching between ASCTB and Azimuth.
    azimuth_colnames <-  colnames(asct_table_derived_from_azimuth)
    if (organ.1=='Brain'){
      azimuth.entire_set_of_cell_types <- asct_table_derived_from_azimuth[grepl("AS/[0-9]$",azimuth_colnames)]
      azimuth.cleaned_set_of_cell_types <- get_cleaned_values_from_df(azimuth.entire_set_of_cell_types)
      n_matching_ct_ontology_ids.4 <- length(intersect(asctb.cleaned_set_of_cell_types, azimuth.cleaned_set_of_cell_types))
    }else{
      azimuth.entire_set_of_ct_ontology_ids <- asct_table_derived_from_azimuth[grepl("AS",azimuth_colnames) & grepl("ID",azimuth_colnames)]
      azimuth.cleaned_set_of_ct_ontology_ids <- get_cleaned_values_from_df(azimuth.entire_set_of_ct_ontology_ids)
      n_matching_ct_ontology_ids.4 <- length(intersect(asctb.cleaned_set_of_ct_ontology_ids, azimuth.cleaned_set_of_ct_ontology_ids))
    }
    
    
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
    return(paste0("Computed statistics for ",config$asctb_name,"."))
  },
  error = function(e){
    cat(paste('\nSomething went wrong while generating the ASCTB organ-level summary stats for:',config$name,'\n'))
    print(e)
  })
}
