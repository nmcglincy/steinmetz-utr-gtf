#
# automatically sorts resulting bed
gtf2bed < hybrid.gtf > hybrid.bed

awk '$5 == "+"' hybrid.bed > hybrid.pos.bed
awk '$5 == "-"' hybrid.bed > hybrid.neg.bed

# bedops --partition hybrid.pos.bed \
# 	| bedmap --count --echo-map - hybrid.pos.bed | head
# looks like id is in the wrong column somehow
# 
# awk -F"|" '{print $1}' hybrid.pos.count.bed > hyb.pos.counts
# 
# bedops --partition hybrid.pos.bed \
# 	| bedmap --count --echo-map - hybrid.pos.bed \
# 	| awk -F"|" '($1 > 2) {print $0}' > weird.ones
# 
# partition and count overlaps
bedops --partition hybrid.pos.bed \
	| bedmap --count --echo - hybrid.pos.bed \
	| awk -F"|" '($1 <2) {print $2}' \
	| bedmap --skip-unmapped --echo-map-range --echo hybrid.pos.bed - \
	| awk -F'["|""\t"]' 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$8,$9,$10,$11,$13}' > hyb-pos-part-simple.bed

# When it got the overlaps it wrote the line from the overlap to file, with the old coordinates...

bedops --partition hybrid.pos.bed \
	| bedmap --count --echo-map - hybrid.pos.bed \
	| awk -F"|" '($1 > 2) {print $2}' > hyb-pos-part-complex.bed

bedops --partition hybrid.neg.bed \
	| bedmap --count --echo - hybrid.neg.bed \
	| awk -F"|" '($1 <2) {print $2}' \
	| bedmap --skip-unmapped --echo-map-range --echo hybrid.pos.bed - \
	| awk -F'["|""\t"]' 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$8,$9,$10,$11,$13}' > hyb-neg-part-simple.bed

bedops --partition hybrid.neg.bed \
	| bedmap --count --echo-map - hybrid.neg.bed \
	| awk -F"|" '($1 > 2) {print $2}' > hyb-neg-part-complex.bed

# I think the best thing might be to reconstruct the mitochondrial genes myself
# Actually, they sorted by location rather than genic identity, so they're all 
# out of order and look weird, so they should be fine as they are actually.

# manually went through complex cases to make *2.bed files

# testing idea
# bedops --partition hybrid.pos.bed \
# 	| bedmap --count --echo-map - hybrid.pos.bed > bar.bed
# 	| awk -F"|" '($1 <2) {print $2}' > foo.bed

# bedmap --skip-unmapped --echo-map-id hybrid.pos.bed hyb-pos-part-simple.bed > foo.bed

# bedmap --skip-unmapped --echo hybrid.pos.bed hyb-pos-part-simple.bed > foo.bed
# 
# putting together the final exon .beds, by strand:
bedops --everything hyb-pos-part-simple.bed \
	mito-pos-ambig.bed \
	hyb-pos-part-complex2.bed > mcglincy.exon.pos.bed

bedops --everything hyb-neg-part-simple.bed \
	hyb-neg-part-complex2.bed > mcglincy.exon.neg.bed
# 
# First I need to extract the CDS lines from our gtf file
# and convert this into a bed.

grep -w "CDS" sac_cer_yassour_utr.gtf | awk -f make-normal-gene-names.awk > sac_cer_yassour_cds.gtf

# wc -l; it's a little short @ 6970, there are some things in there without an ORF?

gtf2bed < sac_cer_yassour_cds.gtf > yassour.cds.bed
awk '$5 == "+"' yassour.cds.bed > yassour.cds.pos.bed
awk '$5 == "-"' yassour.cds.bed > yassour.cds.neg.bed

# ok, I can see one example in the head where removing the exonic overlap had truncated the ORF.
# I think it will be easier of me to think about it in R, but first I will merge and sort everything
# together.

bedops --everything mcglincy.exon.pos.bed \
	mcglincy.exon.neg.bed \
	yassour.cds.bed > needs-to-be-fixed.bed

say -v Milena "I'm finished" 

# there's a bloody formatting error in the gene names 