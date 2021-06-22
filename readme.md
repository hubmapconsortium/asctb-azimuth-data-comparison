# Comparing Azimuth data with ASCT+B reporter data

## Scripts

### ```extract_ct_from_ref.R```

#### Description
This script takes [reference files from Azimuth](https://azimuth.hubmapconsortium.org/references/)
and converts them into csv files that can be viewed in [the ASCT+B reporter](https://hubmapconsortium.github.io/ccf-asct-reporter/).

#### Usage
```
Rscript extract_ct_from_ref.R <organ_name> <hierarchy_column-1> <hierarchy_column-2> <hierarchy_column-3> ...
```

```<organ_name>``` is a file called ```<organ_name>.Rds``` stored in
```azimuth_references``` directory. More organ reference files can be added from
[the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/).

```<hierarchy_column-n>``` are names of columns in the reference Rds files that
define the cell type hierarchy of a cell in the organ. These column names can be
found in the "Annotation Details" section of an organ on [the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/)

The organ's resulting ASCT+B table format file is stored in ```asctb_tables``` directory as
```<organ_name>.csv```. These csv files can be uploaded on google sheets and 
loaded into ASCT+B reporter's playground interface.

#### Example

```
Rscript extract_ct_from_ref.R lungs annotation.l1 annotation.l2
```
