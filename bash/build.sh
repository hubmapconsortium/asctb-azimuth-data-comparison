#! /usr/bin/env bash

mkdir dist
cp data/asctb_tables/*.csv dist
cp data/azimuth_summary_tables/*.csv dist
cp public/* dist
mv dist/index-template.md dist/index.md
# Parse the references[] json from 'organ_data.json' and print each "- organ_name: URL for file_name.csv". Currently hardcoded the string for spleen since we don't have Cell-type info in the custom JSON file on S3.
jq -r '.references[]|"\(.name)|\(.asctb_name)"' data/organ_data.json |\
  awk -F "|" 'BEGIN{OFS=""} {print ($1 ~ /spleen/ ? "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv)" : "- ["$2"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".csv) : [Summary of Cell-Types](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/"$1".celltype_stats.csv)")}'|\
   sort > dist/csv_download_list.md

cat dist/csv_download_list.md >> dist/index.md

cd dist
echo -e "\n## Summary statistics of all datasets processed currently:" >> index.md
ls *all*stat* | awk -F "." 'BEGIN{OFS=""} {print "- [",toupper(substr($1,1,1)),substr($1,2),"](https://hubmapconsortium.github.io/asctb-azimuth-data-comparison/",$1,".",$2,".csv)"}' >> index.md