#!/usr/bin/env bash
# Local (laptop) stages of the v4 analysis, in the order they were run, with the exact arguments.
# The cluster stages (similarity cache, all-pairs decomposition, gene recovery, per-term tests, place-of-care
# contrasts, gene ranking against annotations, variance decomposition) are the sbatch files in this directory
# (allpairs_um, rec_um, terms_um + terms_tasks.txt, setting_um + setting_tasks.txt, vsrest_um, gda_um, r2_um).
#
# DATA holds the inputs that are not in git (see README, "Not in git"):
#   $DATA/hpo/hp.obo, $DATA/hpo/phenotype.hpoa, $DATA/hpo/genes_to_phenotype.txt   (data/get_hpo.sh)
#   $DATA/ontofunc_hp/{term.txt,term2term.txt,graph_path.txt}                     (OntoFUNC + ELK build of hp.owl)
#   $DATA/abs/abstracts.jsonl                                                      (pubmed_abstracts.py output; not shipped)
#   $DATA/hgnc/hgnc_complete_set.txt                                               (HGNC complete set)
#   $DATA/naturalearth/ne_110m_admin_0_countries.shp                               (Natural Earth, public domain)
#   $DATA/storeabs/cands.jsonl                                                     (abstracts of the Store publications, pubmed_abstracts.py)
set -euo pipefail
cd "$(dirname "$0")/.."
DATA=${DATA:-data}
F=results/v4/fine2; S=results/v4/setting2
T=results/termclass.jsonl; G=results/groups_lastauthor.tsv
PY=${PY:-python3}

# 1. Labels: six populations, age and sex, ancestry x place of care
$PY fine_regions.py --packets results/phenopackets --cases results/cases_store.tsv --cases results/cases_pmcp_map2.tsv \
    --cases results/cases_pavs.tsv --cases results/cases_abs.tsv --out-dir results/v4/fine \
    --author-country phenopacket_store_0.1.26_pmid_country.tsv --author-country results/pmcp_pmid_country.tsv \
    --author-country results/pavs_pmid_country.tsv --author-country $DATA/abs/abstracts.jsonl
$PY add_age_sex.py --index results/phenopackets/index.tsv --cases results/v4/fine/cases_store.tsv \
    --cases results/v4/fine/cases_pmcp_map2.tsv --cases results/v4/fine/cases_abs.tsv --out-dir $F
$PY setting_labels.py --packets results/phenopackets --cases $F/cases_abs.tsv --cases $F/cases_pmcp_map2.tsv --out-dir $S

# 2. Tables collected from the cluster outputs (all-pairs and within-vs-rest contrasts, outcome tests, setting
#    contrasts, gene-ranking statistics), then Benjamini-Hochberg q-values per family
$PY collect_tables.py --res results/v4
$PY qvalues.py results/v4 results/v4/qvalues.tsv

# 3. Signs versus symptoms; age and sex interactions; pooled place of care
mkdir -p results/v4/symptom_sign results/v4/age_sex results/v4/setting_pooled
$PY symptom_sign.py --cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv \
    --setting $S/cases_abs.tsv --setting $S/cases_pmcp_map2.tsv --termclass $T --groups $G --perms 2000 \
    --out results/v4/symptom_sign/symptom_sign > results/v4/symptom_sign/log.txt
$PY age_sex_interaction.py --cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv \
    --termclass $T --groups $G --perms 2000 --out results/v4/age_sex/age_sex_interaction > results/v4/age_sex/log.txt
$PY setting_pooled.py --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv --setting $S/cases_abs.tsv \
    --setting $S/cases_pmcp_map2.tsv --termclass $T --groups $G --perms 2000 \
    --out results/v4/setting_pooled/setting_pooled > results/v4/setting_pooled/log.txt

# 4. Organ-system co-reporting networks (all findings, and test-based / bedside split)
mkdir -p results/v4/comorbidity
C="--cases $F/cases_store.tsv --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv --obo $DATA/hpo/hp.obo \
   --termdb $DATA/ontofunc_hp --groups $G --termclass $T --perms 2000"
$PY comorbidity_net.py $C --out results/v4/comorbidity/organ > results/v4/comorbidity/organ.log
$PY comorbidity_net.py $C --split --out results/v4/comorbidity/organ_split > results/v4/comorbidity/organ_split.log

# 5. Country-level covariates (World Bank indicators are downloaded into results/v4/worldbank by the script)
$PY healthspend.py --index results/phenopackets/index.tsv --cases $F/cases_abs.tsv --cases $F/cases_pmcp_map2.tsv \
    --termclass $T --wb-dir results/v4/worldbank --years-abs $DATA/abs/abstracts.jsonl \
    --years-tsv results/pmcp_pmid_country.tsv --out results/v4/worldbank/healthspend

# 6. Gene ranking against annotations: statistics over the cluster chunks (gda_um)
$PY gda_impact.py stats --percase 'results/v4/gda/gda_chunk*.csv' --cases $F/cases_store.tsv \
    --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv --groups $G --perms 2000 --out results/v4/gda/gda \
    > results/v4/gda/stats.log

# 7. Extraction performance against Phenopacket Store curation, and replication by corpus
mkdir -p results/v4/validation
$PY abs_to_cases.py --abstracts $DATA/storeabs/cands.jsonl --hgnc $DATA/hgnc/hgnc_complete_set.txt \
    --out results/v4/validation/cases_storeabs.tsv
$PY validate_against_store.py --store-cases $F/cases_store.tsv --pred-cases results/v4/validation/cases_storeabs.tsv \
    --obo $DATA/hpo/hp.obo --termclass $T --out results/v4/validation/abs_vs_store
$PY validate_against_store.py --store-cases $F/cases_store.tsv --pred-cases $F/cases_pmcp_map2.tsv \
    --obo $DATA/hpo/hp.obo --termclass $T --out results/v4/validation/pmcp_vs_store
$PY corpus_replication.py --fine2 $F --units results/v4/validation/abs_vs_store_units.csv --termclass $T --groups $G \
    --perms 1000 --out results/v4/validation/replication

# 8. Sensitivity: Pakistan and Afghanistan assigned to South Asia
mkdir -p results/v4/sensitivity
$PY sensitivity_pakistan.py --index results/phenopackets/index.tsv --cases $F/cases_store.tsv \
    --cases $F/cases_pmcp_map2.tsv --cases $F/cases_abs.tsv --termclass $T --groups $G --perms 2000 \
    --out results/v4/sensitivity/pakistan > results/v4/sensitivity/pakistan.log

# 9. Variance decomposition of within-gene similarity (cluster chunks from jobs/r2_um.sbatch) and descriptive counts
$PY r2_summary.py --res results/v4
$PY summaries.py --res results/v4

# 10. Figures
FI=results/v4/fine
$PY network_figs.py popnet --print-layout --contrasts results/v4/allpairs/contrasts.tsv --within results/v4/qvalues.tsv \
    --cases $FI/cases_store.tsv --cases $FI/cases_pmcp_map2.tsv --cases $FI/cases_abs.tsv --cases $FI/cases_pavs.tsv \
    --out results/v4/figures/gim_fig4_popnet
$PY network_figs.py shift --res results/v4 --out results/v4/figures/fig_shift
$PY gim_figures.py --res results/v4 --shapes $DATA/naturalearth/ne_110m_admin_0_countries.shp
