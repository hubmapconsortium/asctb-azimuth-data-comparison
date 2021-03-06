#! /usr/bin/env bash

mkdir dist
cp data/asctb_formatted_azimuth_data/*.csv dist

cp data/summary_tables/*.csv dist
cp public/* dist

mkdir dist/staging_asctb_master_data/
cp data/staging_area/*.csv dist/staging_asctb_master_data

mv dist/index-template.md dist/index.md

# Parse the references[] json from 'azimuth_asctb_comparison_config.json' and print each "- organ_name: URL for file_name.csv". Currently hardcoded the string for spleen since we don't have Cell-type info in the custom JSON file on S3.
jq -r '.references[]|"\(.name)|\(.asctb_name)"' data/azimuth_asctb_comparison_config.json |\
  awk -F "|" 'BEGIN{OFS=""} {print ($1 ~ /spleen/ ? "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv)" : "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv) --- [Summary of Cell-Types](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".celltype_stats.csv)")}'| sort >> dist/index.md

echo -e "\n## Summary statistics of all datasets processed currently:" >> dist/index.md

cp data/summary_tables/*.xlsx dist
ls dist | grep summaries | awk -F "." 'BEGIN{OFS=""} {print "- [",toupper(substr($1,1,1)),substr($1,2)," ",$2,"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/",$1,".",$2,".",$3,")"}' | sort -r >> dist/index.md
