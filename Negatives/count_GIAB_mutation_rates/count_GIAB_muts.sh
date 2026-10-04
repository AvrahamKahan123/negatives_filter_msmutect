#!/bin/bash

awk 'BEGIN { OFS="\t"; print "sample", "M_count" }
FNR==1 {
    if (s != "") print s, (found ? m : 0)
    s = FILENAME; sub(/.*\//, "", s); sub(/\.full\.call_counts\.tsv$/, "", s)
    found = 0
}
$2 == "M" && $3 == "(Mutation)" { m = $1; found = 1 }
END { if (s != "") print s, (found ? m : 0) }
' /home/avraham/MaruvkaLab/msmutect_postprocessing/Negatives/count_GIAB_mutation_rates/results/with_fisher_fixed/*.full.call_counts.tsv  > M_counts_GIAB.tsv
