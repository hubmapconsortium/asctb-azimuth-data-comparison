# Comparing Azimuth data with ASCT+B reporter data

## Scripts

### ```R/extract_ct_from_ref.R```

#### Description
This script takes [reference files from Azimuth](https://azimuth.hubmapconsortium.org/references/)
and converts them into csv files that can be viewed in [the ASCT+B reporter](https://hubmapconsortium.github.io/ccf-asct-reporter/).
We also generate summaries for each organ to denote Cell-Type vs Number-of-cells.
This final summary excel sheet available [here](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/) will have:
<ol>
<li>Stats of all organ-datasets processed from Azimuth reference data and ASCT+B master data (Number of CTs, Number of Biomarkers, etc.).</li>
<li>Comparison of Azimuth vs ASCT+B (CTs in Azimuth that are missing in ASCT+B, Biomarkers in Azimuth missing in ASCT+B).</li>
<li>A listing of such CTs and Biomarkers identified in previous point.</li>
</ol>

![Repository documentation diagram](/public/Documentation_Diagrams.jpg?raw=True)

#### Usage
```
Rscript R/extract_ct_from_ref.R
```

This script reads a JSON config file called ```azimuth_asctb_comparison_config.json``` in the
```data/``` directory. Then, based on this config it extracts the data from the
required Rds files and stores them in ASCT+B table format. The structure of this
config file is as follows:

```json
{
  "references": [
    {
      "name": "<organ_name>",
      "url": "<azimuth_reference_data_download_url>",
      "cell_type_columns": [
        "<hierarchy_column-1>",
        "<hierarchy_column-2>",
        ...
        "<hierarchy_column-n>"
      ],
      "asctb_name": "<deployment_display_name>",
	  "asctb_master_url": "<asctb_master_data_download_url>"
    },
    {
      "name": "<organ_name>",
      "url": "<azimuth_reference_data_download_url>",
      "cell_type_columns": [
        "<hierarchy_column-1>",
        "<hierarchy_column-2>",
        ...
        "<hierarchy_column-n>"
      ],
	  "new_cell_type_files":[
		"<csv_file-1>",
		"<csv_file-2>",
		...
		"<csv_file-n>"
      "asctb_name": "<deployment_display_name>",
	  "asctb_master_url": "<asctb_master_data_download_url>"
    },
    {
      ...
    },
    ...
  ] 
}

```

```<azimuth_reference_data_download_url>``` is the URL for downloading the Azimuth 
reference file for ```<organ_name>```. Such download URLs can be found by 
clicking on the Zenodo URLs on the [the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/).
```<organ_name>``` is the name that is associated with its respective download URL.

```<hierarchy_column-n>``` are names of columns in the reference Rds files that
define the cell type hierarchy of a cell in the organ. These column names can be
found in the "Annotation Details" section of an organ on 
[the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/).
Columns higher in the cell type hierarchy should have a lower value of "n".

```<csv_files-n>``` are names of csv-files to be retrieved from the Azimuth Website backend,
that contain the cell-type hierarchy data of a cell in the organ.
These file names match with the annotations mentioned as per the above "hierarchy_col".
File-Names that are higher in the cell type hierarchy should have a lower value of "n".

```<asctb_master_data_download_url>``` is the URL for downloading the v1 Google Sheet 
for each ```<organ_name>```. These download URLs can be found on [the CCF Master Tables page](https://hubmapconsortium.github.io/ccf-asct-reporter/).

```<deployment_display_name>``` is the string that is displayed on the final GitHub 
Pages website for the download link of the corresponding config.

Furthermore, this script also merges ontology ID and ontology label of every cell type
with the cell type data extracted from Azimuth's reference files. This data can be
found in a table on [the Azimuth references webpage](https://azimuth.hubmapconsortium.org/references/)
under the "Annotation Details" section for every available body organ. 

This ontology-data is retrieved by pulling and parsing files from [the Azimuth Website Repository](https://github.com/satijalab/azimuth_website).

The organ's final ASCT+B table format file is stored in ```data/asctb_formatted_azimuth_data``` 
directory as ```<organ_name>.csv```. 

The corresponding summary file is stored in ```data/summary_tables``` 
directory as ```<organ_name>.celltype_stats.csv```.

The overall summary file for ASCT+B master data, and Azimuth reference data ingested 
is stored in ```data/summary_tables``` directory as ```Azimuth_vs_ASCTB.summaries.xlsx```. 
This summary file for all organs will contain the number of matches found for CellTypes and Biomarkers, 
between ASCT+B and Azimuth. It will also identify CellTypes and Biomarkers which are present in Azimuth but not in ASCT+B.
The CellTypes in ASCT+B not present in Azimuth mostly are CellTypes with **similar** names but different CT-IDs, hence causing a mismatch.
Some discrepancies seem to be with the different ways of reporting same data.
Since we can’t rely on exact matches of ASCTB.CTname against Azimuth.CTname, we use the standardized CT-IDs to compare.
The missing-CTs list is = set(Azimuth CT IDs) minus set(ASCT+B CT IDs).

For example-
```
“Classical Dendritic” is present in Azimuth yes, with CT-ID “CL:0000451”.

But the CT with similar name in ASCT+B identified as [missing in Azimuth] is “Dendritic Cell (classical)” has CT-ID “CL:0000990”.

“CL:0000990” doesn’t show up in the Azimuth Kidney annotation files, hence marked as missing-CT.
```

On Ontobee, this [Azimuth CT-ID](http://www.ontobee.org/ontology/CL?iri=http://purl.obolibrary.org/obo/CL_0000451) is at a higher class-hierarchy than the [ASCTB CT-ID](http://www.ontobee.org/ontology/CL?iri=http://purl.obolibrary.org/obo/CL_0000990).


The R-package openxlsx has an issue, while writing HYPERTEXT formula to file.

For now implemented hardcoded formula, so:
1. Double click on the HYPERTEXT formula.

2. Click on another cell. This forces Excel to recognize the HYPERTEXT formula.

3. Click on hyperlink to go to a sheet.

The list o Biomarkers present in Azimuth but not in ASCT+B, has HGNC-IDs retrieved from the API maintained by
the Hugo Gene Nomenclature Committee, with documentation available at [HGNC REST web-service docs](https://www.genenames.org/about/guidelines/).
For identifying this list, we first compare Azimuth_BG_names vs ASCT+B_BG_names. Then the HGNC-ID is retrieved for this list 
of Azimuth_BG_names not present in ASCT+B, for comparison with ASCT+B_BG_IDs.
A local cache-file ```biomarker_name_vs_id_cached.csv``` is maintained in the ```data/``` directory, and is 
overwritten at the end of each run to store HGNC-IDs for the list of Azimuth_BG_names not present in ASCT+B.


## Deployment

All available CSV files can be seen [here](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/).
Through GitHub Workflows every time a commit is pushed to main, the script 
```R/extract_ct_from_ref.R``` is run and resultant CSV files are deployed so that
they can be accessed through GitHub pages. This ensures latest ASCT+B CSV data and summaries
are always available on GitHub pages.
