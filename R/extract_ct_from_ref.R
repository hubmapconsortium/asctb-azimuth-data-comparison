library(Seurat)
library(rjson)


# args <- commandArgs(trailingOnly = TRUE)
#
# body_organ <- args[1]
# cell_hierarchy_cols <- args[2:length(args)]

extract_reference  <- function(body_organ, cell_hierarchy_cols) {

  ref <- readRDS(stringr::str_interp("data/azimuth_references/${body_organ}.Rds"))

  cell_types <- ref@meta.data[cell_hierarchy_cols];
  final_column_names <- sapply(1:length(cell_hierarchy_cols), function (x) return(stringr::str_interp("AS/${x}")))
  names(cell_types) <- final_column_names
  count_col <- stringr::str_interp("AS/${length(cell_hierarchy_cols)}/count")
  cell_types[count_col] <- 1;
  group_by_formula <- as.formula(
    paste(
      paste('`',count_col,'` ~ ',sep = ''),
      paste(sapply(final_column_names,
                   function(col)
                     paste("`",col,"`", sep = '')),
            collapse = " + ")))
  cell_types <- aggregate( group_by_formula, data = cell_types, FUN = sum )

  empty_rows <- matrix(rep(NA, length(cell_hierarchy_cols)*10), 10,length(cell_hierarchy_cols))
  write.table(empty_rows,
              file =  stringr::str_interp("data/asctb_tables/${body_organ}.csv"),
              sep = ',',
              na = "",
              row.names = FALSE,
              col.names = FALSE)
  write.table(cell_types,
              file = stringr::str_interp("data/asctb_tables/${body_organ}.csv"),
              sep = ',',
              append = TRUE,
              row.names = FALSE,
              col.names = TRUE)

}

for (reference in fromJSON(file = 'data/organ_data.json')$references) {
  extract_reference(reference$name, reference$cell_type_columns);
}


