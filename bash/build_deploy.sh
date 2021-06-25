#! /usr/bin/env bash

mkdir deploy
cp data/asctb_tables/*.csv deploy
cp public/* deploy
mv deploy/index-template.md deploy/index.md
jq -r '.references[].name' data/organ_data.json |\
  awk 'BEGIN{OFS=""} {print "- [",$1,"](",$1,".csv)"}' >> deploy/index.md
