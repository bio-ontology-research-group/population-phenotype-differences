set -u
cd /home/leechuck/Documents/papers/pavs-saudi-specific-phenotypes/analysis
P=/home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/statsenv/bin/python
T="--obo /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/hpo/hp.obo --termdb /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/ontofunc_hp --groups results/groups_lastauthor.tsv --termclass results/termclass.jsonl"
term() { nice -n 15 $P population_model.py terms "$@" $T; }
for r in MENA EASIA; do
  for st in total patient-carried; do
    echo "########## [v4 all corpora] $r vs EUR terms ($st)"
    term --cases results/cases_store.tsv --cases results/cases_pmcp_map2.tsv --cases results/cases_abs.tsv --r1 $r --r2 EUR --strata $st --perms 2000 --show 10 --out results/v4/terms_all_${r}_$st
    echo "########## [v4 abstracts only] $r vs EUR terms ($st)"
    term --cases results/cases_abs.tsv --r1 $r --r2 EUR --strata $st --perms 2000 --show 10 --out results/v4/terms_abs_${r}_$st
  done
done
S="--cases results/v4/cases_abs_setting.tsv --cases results/v4/cases_pmcp_map2_setting.tsv"
for pair in "MENAanc_EURcare EURorlocal_EURcare" "MENAanc_MENAcare MENAanc_EURcare" "local_MENAcare MENAanc_MENAcare" "EASIAanc_EURcare EURorlocal_EURcare" "EASIAanc_EASIAcare EASIAanc_EURcare"; do
  set -- $pair
  echo "########## [v4 setting] $1 vs $2 terms (total)"
  term $S --r1 $1 --r2 $2 --strata total --perms 2000 --min-cases 5 --show 10 --out results/v4/terms_setting_$1_vs_$2
done
for f in age_neonate+infant age_child age_adolescent age_adult sex_male sex_female; do
  for r in MENA EASIA; do
    echo "########## [v4 stratum $f] $r vs EUR terms (total)"
    term --cases results/v4/cases_abs_$f.tsv --r1 $r --r2 EUR --strata total --perms 2000 --show 5 --out results/v4/terms_${f}_$r
  done
done
echo ALLDONE
