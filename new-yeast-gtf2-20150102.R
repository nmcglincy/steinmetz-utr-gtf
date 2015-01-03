# MY SECOND ATTEMPT TO INCORPORATE THE STEINMETZ DATA INTO A GTF
# 20150102
# 
# INCLUDE CLASSES: "Covering one intact ORF" and "Intergenic transcripts"
# REORDER SO THAT THE START COORDINATE IS ALWAYS SMALLER THAN THE END COORDINATE
# ONLY INCLUDE THE YPD COUNTS
# 
# 
system("grep -w 'Covering one intact ORF' S2_tcd_mTIFAnno.txt | awk -f reformater.awk > mtifs.txt")
system("grep -w 'Intergenic transcripts' S2_tcd_mTIFAnno.txt | awk -f reformater-igt-tifs.awk >> mtifs.txt")
system("awk -f reorder-loc.awk mtifs.txt > mtifs2.txt")
