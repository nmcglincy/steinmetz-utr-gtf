# !/bin/bash
# 
# reformatting chr names from roman numerals to 0n numbers:
sed -i .tmp -e 's/chrXVI/chr16/g' \
			-e 's/chrXV/chr15/g' \
			-e 's/chrXIV/chr14/g' \
			-e 's/chrXIII/chr13/g' \
			-e 's/chrXII/chr12/g' \
			-e 's/chrXI/chr11/g' \
			-e 's/chrX/chr10/g' \
			-e 's/chrIX/chr09/g' \
			-e 's/chrVIII/chr08/g' \
			-e 's/chrVII/chr07/g' \
			-e 's/chrVI/chr06/g' \
			-e 's/chrV/chr05/g' \
			-e 's/chrIV/chr04/g' \
			-e 's/chrIII/chr03/g' \
			-e 's/chrII/chr02/g' \
			-e 's/chrI/chr01/g' \
			-e 's/chrM/chrmt/g' $1 
# 
# reformatting gene names
awk -F'[\t\" ]' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$11}' $1 > tmp && mv tmp $1
# 
# say -v Milena "Im finished" 