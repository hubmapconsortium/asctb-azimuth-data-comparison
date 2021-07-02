library(Seurat)
library(rjson)

process_reference <- function(organ_config) {
  body_organ <- organ_config$name
  cell_hierarchy_cols <- organ_config$cell_type_columns
  ref_url <- organ_config$url

  file_name <- stringr::str_interp("data/azimuth_references/${body_organ}.Rds")
  if (file.exists(file_name)) {
    ref <- readRDS(file_name)
  } else {
    ref <- readRDS(gzcon(url(ref_url)))
    saveRDS(ref, file_name)
  }

  cell_types <- ref@meta.data[cell_hierarchy_cols];
  final_column_names <- sapply(1:length(cell_hierarchy_cols),
                               function (x)
                                 return(stringr::str_interp("AS/${x}")))
  names(cell_types) <- final_column_names
  count_col <- stringr::str_interp("AS/${length(cell_hierarchy_cols)}/COUNT")
  cell_types[count_col] <- 1;
  group_by_formula <- as.formula(
    paste(
      paste('`',count_col,'` ~ ',sep = ''),
      paste(sapply(final_column_names,
                   function(col)
                     paste("`",col,"`", sep = '')),
            collapse = " + ")))
  cell_types <- aggregate( group_by_formula, data = cell_types, FUN = sum )

  return(cell_types)
}

process_cell_type_data <- function(body_organ, cell_hierarchy_cols) {
  return(sapply(1:length(cell_hierarchy_cols),
                function (i, levels) {
                  ont_data <- read.csv(
                    stringr::str_interp(
                      "data/azimuth_ct_tables/${body_organ}__${levels[i]}.csv"))
                  df <- data.frame(
                    name = ont_data$Label,
                    label = if("OBO.Ontology.ID" %in% names(ont_data))
                              gsub("^\\[(.*?)\\].*", "\\1",ont_data$OBO.Ontology.ID)
                            else rep(NA, nrow(ont_data)),
                    id = if("OBO.Ontology.ID" %in% names(ont_data))
                          gsub(".*http:\\/\\/purl\\.obolibrary\\.org\\/obo\\/(.*?)_(.*?)\\)$",
                               "\\1:\\2",
                               ont_data$OBO.Ontology.ID)
                          else rep(NA, nrow(ont_data)))
                  names(df) <- c(paste("AS/",i, sep = ""),
                                 paste("AS/",i,"/LABEL", sep = ""),
                                 paste("AS/",i,"/ID", sep = ""))
                  return(df)
                  },
                cell_hierarchy_cols,
                simplify = FALSE))
}

write_asctb <- function(body_organ, asctb_table) {
  column_count <- ncol(asctb_table)
  header_rows <- matrix(c(paste(body_organ,
                                "cell types as anatomical structures",
                                "from Azimuth reference data"), rep(NA, column_count-1),
                          rep(NA, column_count),
                          "Author Name(s):", rep(NA, column_count-1),
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
              file =  stringr::str_interp("data/asctb_tables/${body_organ}.csv"),
              sep = ',',
              na = "",
              row.names = FALSE,
              col.names = FALSE)
  write.table(asctb_table,
              file = stringr::str_interp("data/asctb_tables/${body_organ}.csv"),
              sep = ',',
              na = "",
              append = TRUE,
              row.names = FALSE,
              col.names = TRUE)
}

# main loop runs for rach organ in the JSON config
for (reference in fromJSON(file = 'data/organ_data.json')$references) {

  # extract azimuth reference data
  reference_table <- process_reference(reference)
  # extract cell type ontology data provided by Jaison
  ct_ontology_tables <- process_cell_type_data(reference$name, reference$cell_type_columns)
  # map ontology ID and LABELS to reference cell types
  merged_data <- Reduce(merge, ct_ontology_tables, reference_table)
  # reorder columns
  column_order <- c(
    sapply(1:length(reference$cell_type_columns),
           function (n)
             c(paste("AS/",n,sep=""),
               paste("AS/",n,"/LABEL",sep=""),
               paste("AS/",n,"/ID",sep=""))),
    paste("AS/",length(reference$cell_type_columns),"/COUNT",sep=""))
  # generate final CSV file
  write_asctb(reference$name, merged_data[,column_order])
}


