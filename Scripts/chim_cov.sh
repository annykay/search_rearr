#14-41
grep u[um] SourceData/results_14_41/gmapped_parsed_sorted_chunks/s4114_1.lane1.14-41.0.pairsam  | grep -v \# | cut -f 2,6,11,13 > DerivedData/s4114_1.14-41.u1.not.pairsam 
grep u[um] SourceData/results_14_41/gmapped_parsed_sorted_chunks/s4114_2.lane1.14-41.0.pairsam  | grep -v \# | cut -f 2,6,11,13 > DerivedData/s4114_2.14-41.u1.not.pairsam 
grep u[um] SourceData/results_14_41/gmapped_parsed_sorted_chunks/X44.lane1.14-41.0.pairsam  | grep -v \# | cut -f 2,6,11,13 > DerivedData/X44.14-41.not.pairsam 

grep [um]u SourceData/results_14_41/gmapped_parsed_sorted_chunks/s4114_1.lane1.14-41.0.pairsam  | grep -v \# | cut -f 4,7,12,14 > DerivedData/s4114_1.14-41.u2.not.pairsam 
grep [um]u SourceData/results_14_41/gmapped_parsed_sorted_chunks/s4114_2.lane1.14-41.0.pairsam  | grep -v \# | cut -f 4,7,12,14 > DerivedData/s4114_2.14-41.u2.not.pairsam 
grep [um]u SourceData/results_14_41/gmapped_parsed_sorted_chunks/X44.lane1.14-41.0.pairsam  | grep -v \# | cut -f 4,7,12,14 > DerivedData/X44.14-41.u2.not.pairsam 

cat DerivedData/s4114_1.14-41.u1.not.pairsam DerivedData/s4114_2.14-41.u1.not.pairsam  DerivedData/s4114_1.14-41.u2.not.pairsam DerivedData/s4114_2.14-41.u2.not.pairsam > DerivedData/chimers/chimeric3parts.bed

awk '{print $1"\t"$4"\t"$4}' DerivedData/chimers/chimeric3parts.bed > DerivedData/chimers/chimeric3.bed
bedtools coverage -a DerivedData/chimers/chimeric3.bed -b SourceData/results_14_41/MBoI25.bed > DerivedData/coverages/chimers3_inters.cov
awk '{if ($4==0) print $1"\t"$2"\t"$3}' DerivedData/coverages/chimers3_inters.cov > DerivedData/coverages/chimer3_not_restriction.bed
wc -l DerivedData/coverages/chimer3_not_restriction.bed
bedtools makewindows -g SourceData/indexes/14-41_scaffolds.chrom.sizes -w 20 > DerivedData/coverages/windows.bed

bedtools coverage -a DerivedData/coverages/windows.bed -b DerivedData/coverages/chimer3_not_restriction.bed > DerivedData/coverages/chimers3_not_restr_binarized.cov

bedtools makewindows -g SourceData/indexes/14-41_scaffolds.chrom.sizes -w 40 > DerivedData/coverages/windows40.bed
bedtools coverage -a DerivedData/coverages/windows40.bed -b DerivedData/coverages/chimer3_not_restriction.bed > DerivedData/coverages/chimers3_not_restr_binarized40.cov


cat DerivedData/X44.14-41.not.pairsam DerivedData/X44.14-41.u2.not.pairsam > DerivedData/chimers/alienChimeric3Parts.bed
awk '{print $1"\t"$4"\t"$4}' DerivedData/chimers/alienChimeric3Parts.bed > DerivedData/chimers/alienChimeric3.bed
bedtools coverage -a DerivedData/chimers/alienChimeric3.bed -b SourceData/results_14_41/MBoI25.bed > DerivedData/coverages/alienChimers3_inters.cov
awk '{if ($4==0) print $1"\t"$2"\t"$3}' DerivedData/coverages/alienChimers3_inters.cov > DerivedData/coverages/alienChimer3_not_restriction.bed
bedtools makewindows -g SourceData/indexes/14-41_scaffolds.chrom.sizes -w 88 > DerivedData/coverages/aWindows.bed

bedtools coverage -a DerivedData/coverages/aWindows.bed -b DerivedData/coverages/alienChimer3_not_restriction.bed > DerivedData/coverages/alienChimers3_not_restr_binarized.cov

bedtools coverage -a DerivedData/coverages/windows40.bed -b DerivedData/coverages/alienChimer3_not_restriction.bed > DerivedData/coverages/alienChimers3_not_restr_binarized40.cov

rm DerivedData/*.not.pairsam 
