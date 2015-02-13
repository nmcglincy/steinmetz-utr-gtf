#
# automatically sorts resulting bed
gtf2bed < hybrid.gtf > hybrid.bed

awk '$5 == "+"' hybrid.bed > hybrid.pos.bed
awk '$5 == "-"' hybrid.bed > hybrid.neg.bed

