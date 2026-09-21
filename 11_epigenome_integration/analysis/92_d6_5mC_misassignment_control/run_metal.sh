#!/bin/zsh
# GPU (Metal) basecall of the two validation files with the manuscript's model. Uses the Apple GPU, little CPU. ~5-10 min.
cd "/Users/okaban/.claude-science/orgs/0710cf94-f852-46cd-9bbe-c2e0c51c8927/workspaces/49ff2ed5-8535-42f7-9e3e-6c98b582b12e/d6_control"
D=/Applications/MinKNOW.app/Contents/Resources/dorado
for s in 5mC_rep1 control_rep1; do
  $D/bin/dorado basecaller $D/data/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 $s.pod5 --modified-bases-models $D/data/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_4mC_5mC@v1 --reference all_5mers.fa --device metal > $s.4mC_5mC.bam
done
echo METAL_DONE > metal.flag
