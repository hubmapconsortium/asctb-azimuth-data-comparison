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
      "url": "<reference_data_download_url>",
      "cell_type_columns": [
        "<hierarchy_column-1>",
        "<hierarchy_column-2>",
        ...
        "<hierarchy_column-n>"
      ]
    },
    {
      "name": "<organ_name>",
      "url": "<reference_data_download_url>",
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

```<reference_data_download_url>``` is the URL for downloading the Azimuth 
reference file for ```<organ_name>```. Such download URLs can be found by 
clicking on the Zenodo URLs on the [the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/).
```<organ_name>``` is the name that is associated with its respective download URL.

```<hierarchy_column-n>``` are names of columns in the reference Rds files that
define the cell type hierarchy of a cell in the organ. These column names can be
found in the "Annotation Details" section of an organ on 
[the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/).
Columns higher in the cell type hierarchy should have a lower value of "n".

Furthermore, this script also merges ontology ID and ontology label of every cell type
with the cell type data extracted from Azimuth's reference files. This data can be 
found in a table on [the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/)
under the section for every available body organ. This ontology data is stored in the
```data/azimuth_ct_tables``` directory where each file is named in the format
```<organ_name>__<hierarchy_column-n>.csv```. Unfortunately, at the time of writing 
this document, these CSV files aren't hosted anywhere to be readily downloaded. These
files were provided to us by courtesy of Jaison from the Azimuth team.

The organ's final ASCT+B table format file is stored in ```data/asctb_tables``` 
directory as ```<organ_name>.csv```. 

## Deployment

The available CSV files can be seen [here](https://darshalshetty.github.io/asctb-azimuth-data-comparison/).
Through GitHub Workflows every time a commit is pushed to main, the script 
```R/extract_ct_from_ref.R``` is run and resultant CSV files are deployed so that
they can be accessed through GitHub pages. Thus ensuring latest ASCT+B CSV data
is always available on GitHub pages.
