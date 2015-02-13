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
	| awk -F"|" '($1 <2) {print $2}' > hyb-pos-part-simple.bed

bedops --partition hybrid.pos.bed \
	| bedmap --count --echo-map - hybrid.pos.bed \
	| awk -F"|" '($1 > 2) {print $2}' > hyb-pos-part-complex.bed

bedops --partition hybrid.neg.bed \
	| bedmap --count --echo - hybrid.neg.bed \
	| awk -F"|" '($1 <2) {print $2}' > hyb-neg-part-simple.bed

bedops --partition hybrid.neg.bed \
	| bedmap --count --echo-map - hybrid.neg.bed \
	| awk -F"|" '($1 > 2) {print $2}' > hyb-neg-part-complex.bed

say -v Milena "I'm finished"