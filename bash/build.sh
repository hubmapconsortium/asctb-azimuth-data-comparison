#! /usr/bin/env bash

mkdir dist
cp data/asctb_tables/*.csv dist
cp public/* dist
mv dist/index-template.md dist/index.md
jq -r '.references[].name' data/organ_data.json |\
  awk 'BEGIN{OFS=""} {print "- [",$1,"](",$1,".csv)"}' >> dist/index.md
