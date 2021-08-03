#! /usr/bin/env bash

mkdir dist
cp data/asctb_tables/*.csv dist
cp public/* dist
mv dist/index-template.md dist/index.md
jq -r '.references[]|"\(.name)|\(.asctb_name)"' data/organ_data.json |\
  awk -F "|" 'BEGIN{OFS=""} {print "- [",$2,"](",$1,".csv)"}' >> dist/index.md
