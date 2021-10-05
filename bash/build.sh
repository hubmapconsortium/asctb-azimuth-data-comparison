#! /usr/bin/env bash

mkdir dist
cp data/asctb_tables/*.csv dist
cp data/summary_tables/*.csv dist
cp public/* dist
mv dist/index-template.md dist/index.md

# Parse the references[] json from 'azimuth_asctb_comparison_config.json' and print each "- organ_name: URL for file_name.csv". Currently hardcoded the string for spleen since we don't have Cell-type info in the custom JSON file on S3.

# todo: Uncomment this update after we finalize the logic for cell-type vs count summaries. Comment out the json parsing code on next line when updating.
#jq -r '.references[]|"\(.name)|\(.asctb_name)"' data/azimuth_asctb_comparison_config.json |\
#  awk -F "|" 'BEGIN{OFS=""} {print ($1 ~ /spleen/ ? "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv)" : "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv) : [Summary of Cell-Types](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".celltype_stats.csv)")}'| sort > dist/csv_download_list.md

jq -r '.references[]|"\(.name)|\(.asctb_name)"' data/azimuth_asctb_comparison_config.json |\
  awk -F "|" 'BEGIN{OFS=""} {print ($1 ~ /spleen/ ? "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv)" : "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv)")}'|\
   sort > dist/csv_download_list.md



cat dist/csv_download_list.md >> dist/index.md

cd dist

echo -e "\n## Summary statistics of all datasets processed currently:" >> index.md

ls *All*stat* | awk -F "." 'BEGIN{OFS=""} {print "- [",toupper(substr($1,1,1)),substr($1,2)," ",$2,"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/",$1,".",$2,".",$3,".csv)"}' | sort -r >> index.md
