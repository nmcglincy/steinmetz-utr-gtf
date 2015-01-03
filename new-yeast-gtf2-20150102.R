# MY SECOND ATTEMPT TO INCORPORATE THE STEINMETZ DATA INTO A GTF
# 20150102
# 
# INCLUDE CLASSES: "Covering one intact ORF" and "Intergenic transcripts"
# REORDER SO THAT THE START COORDINATE IS ALWAYS SMALLER THAN THE END COORDINATE
# ONLY INCLUDE THE YPD COUNTS
# 
system("grep -w 'Covering one intact ORF' S2_tcd_mTIFAnno.txt | awk -f reformater.awk > mtifs.txt")
system("grep -w 'Intergenic transcripts' S2_tcd_mTIFAnno.txt | awk -f reformater-igt-tifs.awk >> mtifs.txt")
system("awk -f reorder-loc.awk mtifs.txt > mtifs2.txt")
#
# READ IN THE RESULTING FILE
mtifs = read.delim("mtifs2.txt",
                   header = FALSE)
# head(mtifs)
colnames(mtifs) = c("chr", "strand", "start", "end", "ypd.counts", "class", "gene.id")
# QUICK SANITY CHECK
with(mtifs, hist(log10(ypd.counts)))
# LOOKS LIKE WHAT I SAW EARLIER
