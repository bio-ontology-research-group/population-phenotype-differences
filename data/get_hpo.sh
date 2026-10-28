#!/bin/bash
# HPO release used for every analysis (hp/releases/2026-01-08) into data/hpo/.
set -euo pipefail
R=https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-01-08
mkdir -p "$(dirname "$0")/hpo"; cd "$(dirname "$0")/hpo"
for f in hp.obo hp.owl phenotype.hpoa genes_to_phenotype.txt phenotype_to_genes.txt; do curl -fsSLO "$R/$f"; done
grep -m1 data-version hp.obo
