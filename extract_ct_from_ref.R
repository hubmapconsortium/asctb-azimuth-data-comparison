#library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

body_organ <- args[1]
cell_hierarchy_cols <- args[2:length(args)]

ref <- readRDS(stringr::str_interp("azimuth_references/${body_organ}.Rds"))

cell_types <- unique(ref@meta.data[cell_hierarchy_cols])
# sorting not necessary
# cell_types <- cell_types[order(cell_types$annotation.l1, cell_types$annotation.l2),]
names(cell_types) <- sapply(1:length(cell_hierarchy_cols), function (x) return(stringr::str_interp("AS/${x+1}")))
cell_types$`AS/1` <- as.factor(rep(body_organ , nrow(cell_types)))

empty_rows <- matrix(rep(NA, length(cell_hierarchy_cols)*10), 10,length(cell_hierarchy_cols))
write.table(empty_rows, 
            file =  stringr::str_interp("asctb_tables/${body_organ}.csv"),
            sep = ',',
            na = "",
            row.names = FALSE, 
            col.names = FALSE)
write.table(cell_types[,c(length(cell_hierarchy_cols)+1,
                          1:length(cell_hierarchy_cols))], # reorder columns
            file = stringr::str_interp("asctb_tables/${body_organ}.csv"),
            sep = ',',
            append = TRUE,
            row.names = FALSE, 
            col.names = TRUE)

