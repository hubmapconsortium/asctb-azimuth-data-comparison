# Comparing Azimuth data with ASCT+B reporter data

## Scripts

### ```R/extract_ct_from_ref.R```

#### Description
This script takes [reference files from Azimuth](https://azimuth.hubmapconsortium.org/references/)
and converts them into csv files that can be viewed in [the ASCT+B reporter](https://hubmapconsortium.github.io/ccf-asct-reporter/).

#### Usage
```
Rscript R/extract_ct_from_ref.R
```

This script reads a JSON config file called ```organ_data.json``` in the
```data``` directory. Then, based on this config it extracts the data from the
required Rds files and stores them in ASCT+B table format. The structure of this
config file is as follows:

```json
{
  "references": [
    {
      "name": "<organ_name>",
      "cell_type_columns": [
        "<hierarchy_column-1>",
        "<hierarchy_column-2>",
        ...
        "<hierarchy_column-n>"
      ]
    },
    {
      "name": "<organ_name>",
      "cell_type_columns": [
        "<hierarchy_column-1>",
        "<hierarchy_column-2>",
        ...
        "<hierarchy_column-n>"
      ]
    },
    {
      ...
    },
    ...
  ] 
}

```


```<organ_name>``` is a file called ```<organ_name>.Rds``` stored in
```data/azimuth_references``` directory. More organ reference files can be added
from [the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/).

```<hierarchy_column-n>``` are names of columns in the reference Rds files that
define the cell type hierarchy of a cell in the organ. These column names can be
found in the "Annotation Details" section of an organ on 
[the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/).
Columns higher in the cell type hierarchy should have a lower value of "n".

The organ's resulting ASCT+B table format file is stored in ```data/asctb_tables``` 
directory as ```<organ_name>.csv```. These csv files can be uploaded on google 
sheets and loaded into ASCT+B reporter's playground interface.

## Deployment

The available CSV files can be seen [here](https://darshalshetty.github.io/asctb-azimuth-data-comparison/).
Through GitHub Workflows every time a commit is pushed to main, the script 
```R/extract_ct_from_ref.R``` is run and resultant CSV files are deployed so that
they can be accessed through GitHub pages. Thus ensuring latest ASCT+B CSV data
is always available on GitHub pages.
