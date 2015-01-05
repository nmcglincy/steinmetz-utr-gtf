# MY SECOND ATTEMPT TO INCORPORATE THE STEINMETZ DATA INTO A GTF
# 20150102
# 
# INCLUDE CLASSES: "Covering one intact ORF" and "Intergenic transcripts"
# REORDER SO THAT THE START COORDINATE IS ALWAYS SMALLER THAN THE END COORDINATE
# ONLY INCLUDE THE YPD COUNTS
#
# ONLY REQUIRED THE FIRST TIME YOU RUN THE SCRIPT
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
# with(mtifs, hist(log10(ypd.counts)))
# LOOKS LIKE WHAT I SAW EARLIER
# 
# I THINK I'M GOING TO USE DIFFERENT STRATEGIES FOR THE TWO DIFFERENT MTIF CATEGORIES, SO I'LL 
# SUBSET BY THAT
coio.mtifs = subset(mtifs, class == "Covering_one_intact_ORF")
igt.mtifs = subset(mtifs, class != "Covering_one_intact_ORF")
 
# head(coio.mtifs)
# summary(coio.mtifs)
# head(igt.mtifs)
# summary(igt.mtifs)
# EVERYTHING LOOKS AS I EXPECTED

# INTERGENIC TRANSCRIPTS
# 
# FOR THE INTERGENIC TRANSCRIPTS I THINK IT'S ENOUGH, IN THE FIRST INSTANCE, TO CALCULATE THE UNION
# OF THE IDENTIFIED FRAGMENTS WITHOUT CONSIDERING THE NUMBER OF READS.

# STEP 1. CONVERT DATAFRAME INTO A GR RANGES OBJECT
library("GenomicRanges")
?GenomicRanges
igt.mtifs.gr = makeGRangesFromDataFrame(igt.mtifs,
                                        keep.extra.columns = TRUE,
                                        ignore.strande = FALSE,
                                        seqnames.field = c("chr"),
                                        start.field = c("start"),
                                        end.field = c("end"),
                                        strand.fiels = c("strand"))







