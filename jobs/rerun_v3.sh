set -u
cd /home/leechuck/Documents/papers/pavs-saudi-specific-phenotypes/analysis
P=/home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/statsenv/bin/python
C="--cases results/cases_store.tsv --cases results/cases_pmcp_map2.tsv --cases results/cases_pavs.tsv"
for r in MENA EASIA; do
  echo "########## [v3 data, new mapping] $r vs EUR: decomposition"
  nice -n 15 python3 population_model.py pairs $C --hpo-dir /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/hpo --cache /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/simcache_merged_v3.npy --groups results/groups_lastauthor.tsv --r1 $r --r2 EUR --perms 1000 --boot 200 --out results/popmodel_v3_$r 2>&1 | grep -v "sim rows\|^cases"
  for st in total patient-carried; do
    echo "########## [v3, new mapping] $r vs EUR: terms ($st)"
    nice -n 15 $P population_model.py terms --cases results/cases_store.tsv --cases results/cases_pmcp_map2.tsv --obo /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/hpo/hp.obo --termdb /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/ontofunc_hp --groups results/groups_lastauthor.tsv --termclass results/termclass.jsonl --r1 $r --r2 EUR --strata $st --perms 2000 --show 10 --out results/popmodel_v3_${r}_$st
  done
  echo "########## [v3, OLD mapping, for comparison] $r vs EUR: terms (total)"
  nice -n 15 $P population_model.py terms --cases results/cases_store.tsv --cases results/cases_pmcp.tsv --obo /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/hpo/hp.obo --termdb /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/ontofunc_hp --groups results/groups_lastauthor.tsv --termclass results/termclass.jsonl --r1 $r --r2 EUR --strata total --perms 2000 --show 0 --out /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/oldmap_${r} | grep -E "cases|terms|global|share|physio|laboratory"
done
LIT=store,pmcp,pavs_lit
RC="$C --hpo-dir /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/hpo --cache /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/simcache_merged_v3.npy --authors results/all_authors.tsv --groups results/groups_lastauthor.tsv --perms 1000"
for spec in "MENA:lit:$LIT" "MENA:clinSA:pavs_clin" "EASIA:lit:$LIT"; do
  IFS=: read r n q <<< "$spec"
  echo "########## [v3, new mapping] retrieval $r vs EUR, queries=$q"
  nice -n 15 python3 retrieval_perm.py $RC --r1 $r --r2 EUR --query-corpus $q --lib-corpus $LIT --out results/popretr_v3_${r}_$n 2>&1 | grep -v "sim rows\|^cases"
done
echo ALLDONE
