library(Seurat)
library(rjson)
library(httr)
source('R/extract_ct_from_json.R')

process_reference <- function(organ_config) {
  body_organ <- organ_config$name
  cell_hierarchy_cols <- organ_config$cell_type_columns
  ref_url <- organ_config$url

  # load input reference files locally if present, else download them from URL in config
  file_name <- paste0("data/azimuth_references/", body_organ, ".Rds")
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
}

process_cell_type_data <- function(body_organ, cell_hierarchy_cols) {
  return(sapply(1:length(cell_hierarchy_cols),
                function (i, levels) {
                  ct_file <- paste0("data/azimuth_ct_tables/", body_organ, "__",
                                    levels[i], ".csv")

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

                  # rename columns according to ASCT+B table format
                  names(df) <- c(paste0("AS/",i),
                                 paste0("AS/",i,"/LABEL"),
                                 paste0("AS/",i,"/ID"))
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
              file =  paste0("data/asctb_tables/", body_organ, ".csv"),
              sep = ',',
              na = "",
              row.names = FALSE,
              col.names = FALSE)
  write.table(asctb_table,
              file = paste0("data/asctb_tables/", body_organ, ".csv"),
              sep = ',',
              na = "",
              append = TRUE,
              row.names = FALSE,
              col.names = TRUE)
}

process_config <- function(config) {
  # extract azimuth reference data
  reference_table <- process_reference(config)
  # extract cell type ontology data provided by Jaison
  ct_ontology_tables <- process_cell_type_data(config$name, config$cell_type_columns)
  # map ontology ID and LABELS to reference cell types
  merged_data <- Reduce(merge, ct_ontology_tables, reference_table)
  # reorder columns
  column_order <- c(
    sapply(1:length(config$cell_type_columns),
           function (n)
             c(paste0("AS/",n),
               paste0("AS/",n,"/LABEL"),
               paste0("AS/",n,"/ID"))),
    paste0("AS/",length(config$cell_type_columns),"/COUNT"))
  # generate final CSV file
  return(merged_data[,column_order])
}

# main loop runs for each organ in the JSON config
for (config in fromJSON(file = 'data/organ_data.json')$references) {
  asct_table <- switch(
    config$mode %||% '',
    'nested-json' = { extract_ct_from_json(config$url) },
    process_config(config)
  )
  write_asctb(config$name, asct_table)
}


