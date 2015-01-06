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
# ?GenomicRanges
igt.mtifs.gr = makeGRangesFromDataFrame(igt.mtifs,
                                        keep.extra.columns = TRUE,
                                        ignore.strand = FALSE,
                                        seqnames.field = c("chr"),
                                        start.field = c("start"),
                                        end.field = c("end"),
                                        strand.field = c("strand"))
# igt.mtifs.gr
# STEP 2. SUBSUME OVERLAPPING RANGES.
# ?reduce
# 
# NOT SURE ABOUT THE with.revmap OR min.gapwidth OPTIONS.
# WHAT IS A GOOD BIOLOGICALLY RELEVANT VALUE OF MIN.GAPWIDTH
igt.mtifs.gr.rd = reduce(igt.mtifs.gr,
                         drop.empty.ranges = FALSE,
                         min.gapwidth = 0,
                         with.revmap = TRUE)
# igt.mtifs.gr.rd
# igt.mtifs.gr[52]
# I THINK IT WOULD BE COOL TO LOOK AT THE LAST TWO OBJECTS IN IGV AS A BIT OF A SANITY CHECK.
# 
# ALSO, IF I COULD LOOK AT THE DISTRIBUTION OF INTER-RANGE DISTANCES, THIS WOULD GUIDE MY CHOICE
# OF min.gapwidth.
# mtif.gaps = gaps(igt.mtifs.gr.rd)
# mtif.gaps
# mtif.gaps.wth = width(mtif.gaps)
# hist(mtif.gaps.wth)
# TOO BROAD A RANGE TO BE INFORMATIVE ABOUT THE SMALLEST GAPS
# mtif.gaps.wth.l100 = mtif.gaps.wth[which(mtif.gaps.wth <= 100)]
# hist(mtif.gaps.wth.l100)
# plot(mtif.gaps.wth.l100)
# I THINK I NEED GGPLOT2 FOR THIS
# library("ggplot2")
# mtif.gaps.wth.l100.df = as.data.frame(mtif.gaps.wth.l100)
# mtif.gaps.wth.l100.df
# ggplot(mtif.gaps.wth.l100.df, aes(x = mtif.gaps.wth.l100)) +
#   geom_dotplot(binwidth = 1, stackdir = "center", pch = 21) +
#   ylab(NULL) +
#   xlab("Inter mTIF-assembly gap") +
#   scale_x_continuous(breaks = seq(0, 100, by = 10)) +
#   theme(panel.border = element_rect(fill = NA, colour = "black"),
#         axis.title.x = element_text(vjust = 0, size = 16),
#         axis.title.y = element_text(vjust = 1, size = 14),
#         axis.text.x = element_text(size=14, vjust = 0.5),
#         axis.text.y  = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.title = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12))
# ggsave("mtif-gaps-lt100.png", dpi = 400)
# SHOULD PLOT THE SAME FOR THE WHOLE VECTOR FOR COMPARISON
# mtif.gaps.wth.df = as.data.frame(mtif.gaps.wth)
# ggplot(mtif.gaps.wth.df, aes(x = log2(mtif.gaps.wth))) +
#   geom_histogram() +
#   geom_vline(xintercept = log2(100), colour = "red") +
#   xlab("Log2 (Inter mTIF-assembly gap)") +
#   ylab("Count") +
#   theme(panel.border = element_rect(fill = NA, colour = "black"),
#         axis.title.x = element_text(vjust = 0, size = 16),
#         axis.title.y = element_text(vjust = 1, size = 16),
#         axis.text.x = element_text(size=14, vjust = 0.5),
#         axis.text.y = element_text(size=14, vjust = 0.5),
#         plot.title = element_text(size = 16),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12))
# ggsave("mtif-gaps.png", dpi = 400)
# 
# THE SMALLEST GAPS ARE
# head(sort(mtif.gaps.wth), n = 10)
#  [1]  1  5  5  7  9  9 10 16 18 19
# WHICH GAP
# which.min(mtif.gaps.wth)
# mtif.gaps[157]
# 
# EFFECT OF VARYING MIN.GAPWIDTH
# length(reduce(igt.mtifs.gr, drop.empty.ranges = FALSE, min.gapwidth = 0, with.revmap = TRUE))
# [1] 2998
# length(reduce(igt.mtifs.gr, drop.empty.ranges = FALSE, min.gapwidth = 1, with.revmap = TRUE))
# [1] 2993
# SO THERE ARE 5 MTIF-ASSEMBLIES THAT ARE DIRECTLY ADJACENT TO EACH OTHER
# 
# DO THE SET OPERATION ON THE DATAFRAMES USING DPLYR
# gap0 = reduce(igt.mtifs.gr, drop.empty.ranges = FALSE, min.gapwidth = 0, with.revmap = FALSE)
# gap0
# gap1 = reduce(igt.mtifs.gr, drop.empty.ranges = FALSE, min.gapwidth = 1, with.revmap = FALSE)
# gap1
# gap0.df = as.data.frame(gap0)
# gap0.df
# gap1.df = as.data.frame(gap1)
# gap1.df
# library(dplyr)
# 
# WHAT IS IN GAP0.DF THAT IS NOT IN GAP1.DF
# dplyr::setdiff(gap0.df,gap1.df)
#   seqnames  start    end width strand
#   1        11 259577 259705   129      +
#   2        11 259706 261084  1379      +
#   3        16 911828 912080   253      +
#   4        16 912081 912564   484      +
#   5         2 538924 539552   629      -
#   6         2 539553 540666  1114      -
#   7         5 101296 101504   209      +
#   8         5 101505 101912   408      +
#   9         7 511009 512949  1941      -
#   10        7 512950 514367  1418      -
# 
# WHAT IS IN GAP1.DF THAT IS NOT IN GAP0.DF
# dplyr::setdiff(gap1.df,gap0.df)
#   seqnames  start    end width strand
#   1       11 259577 261084  1508      +
#   2       16 911828 912564   737      +
#   3        2 538924 540666  1743      -
#   4        5 101296 101912   617      +
#   5        7 511009 514367  3359      -
# 
# AFTER THAT, I THINK IT'S THE BEST THING TO TO STIPULATE min.gapwidth = 0

# EACH IGT GETS A UNIQUE ID, WHICH IS IGT_ FOLLOWED BY THE CHR, THE A _, FOLLOWED BY THE ROW NUMBER
# OF THE TABLE. ALSO ADDED THE NUMBER OF YPD READS FOR EACH IGT, MIGHT BE USEFULL TO PRIORITISE
# THINGS LATER.
igt.mtifs.gr.rd.df = as.data.frame(igt.mtifs.gr.rd)
igt.mtif.asm = data.frame(igt.mtifs.gr.rd.df[,1:5],
                          gene.id = sprintf("igt_%s_%d", 
                                            igt.mtifs.gr.rd.df$seqnames, 
                                            1:nrow(igt.mtifs.gr.rd.df)),
                          ypd.counts = unlist(lapply(igt.mtifs.gr.rd.df$revmap, 
                                                     function(x) {sum(igt.mtifs.gr[x]$ypd.counts)})))
igt.mtif.asm$seqnames = sprintf("chr%d", igt.mtif.asm$seqnames)
# head(igt.mtif.asm)
# tail(igt.mtif.asm)
#  - GRAPH MIGHT BE NICE - WIDTH AND NO.READS DISTRIBUTION
library(ggplot2)
ggplot(igt.mtif.asm, aes(x = width)) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(0, 6000, by = 1000)) +  
  ylab("Counts") +
  xlab("IGT assembly width") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y = element_text(size=14, vjust = 0.5),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
ggsave("mtif-asm-wdth.png", dpi = 400)

ggplot(igt.mtif.asm, aes(x = log10(ypd.counts))) +
  geom_histogram() +
  ylab("Counts") +
  xlab("log10(reads in YPD)") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 16),
        axis.title.y = element_text(vjust = 1, size = 16),
        axis.text.x = element_text(size=14, vjust = 0.5),
        axis.text.y = element_text(size=14, vjust = 0.5),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
ggsave("mtif-asm-ypdReads.png", dpi = 400)

# WRTIE OUT FOR POSTERITY
write.csv(igt.mtif.asm,
          file = "igt-mtifs-asm.csv",
          row.names = FALSE,
          quote = FALSE)

# 




