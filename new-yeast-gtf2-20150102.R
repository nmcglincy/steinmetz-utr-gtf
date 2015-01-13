# MY SECOND ATTEMPT TO INCORPORATE THE STEINMETZ DATA INTO A GTF
# 20150102
# 
# INCLUDE CLASSES: "Covering one intact ORF" and "Intergenic transcripts"
# REORDER SO THAT THE START COORDINATE IS ALWAYS SMALLER THAN THE END COORDINATE
# ONLY INCLUDE THE YPD COUNTS
#
# ONLY REQUIRED THE FIRST TIME YOU RUN THE SCRIPT
system("grep -w 'Covering one intact ORF' S2_tcd_mTIFAnno.txt | grep -wv 'NA' | awk -f reformater.awk > mtifs.txt")
system("grep -w 'Intergenic transcripts' S2_tcd_mTIFAnno.txt | awk -f reformater-igt-tifs.awk >> mtifs.txt")
system("awk -f reorder-loc.awk mtifs.txt > mtifs2.txt")
#
# READ IN THE RESULTING FILE
mtifs = read.delim("mtifs2.txt",
                   header = FALSE)
colnames(mtifs) = c("chr", "strand", "start", "end", "ypd.counts", "class", "gene.id")
# 
# I THINK I'M GOING TO USE DIFFERENT STRATEGIES FOR THE TWO DIFFERENT MTIF CATEGORIES, SO I'LL 
# SUBSET BY THAT
coio.mtifs = subset(mtifs, class == "Covering_one_intact_ORF")
igt.mtifs = subset(mtifs, class != "Covering_one_intact_ORF")
# 
# INTERGENIC TRANSCRIPTS
# 
# FOR THE INTERGENIC TRANSCRIPTS I THINK IT'S ENOUGH, IN THE FIRST INSTANCE, TO CALCULATE THE UNION
# OF THE IDENTIFIED FRAGMENTS WITHOUT CONSIDERING THE NUMBER OF READS.
# 
# STEP 1. CONVERT DATAFRAME INTO A GR RANGES OBJECT
library("GenomicRanges")
igt.mtifs.gr = makeGRangesFromDataFrame(igt.mtifs,
                                        keep.extra.columns = TRUE,
                                        ignore.strand = FALSE,
                                        seqnames.field = c("chr"),
                                        start.field = c("start"),
                                        end.field = c("end"),
                                        strand.field = c("strand"))
# 
# STEP 2. SUBSUME OVERLAPPING RANGES.
# 
# NOT SURE ABOUT THE with.revmap OR min.gapwidth OPTIONS.
# WHAT IS A GOOD BIOLOGICALLY RELEVANT VALUE OF MIN.GAPWIDTH
igt.mtifs.gr.rd = reduce(igt.mtifs.gr,
                         drop.empty.ranges = FALSE,
                         min.gapwidth = 0,
                         with.revmap = TRUE)
# 
# 
# ALSO, IF I COULD LOOK AT THE DISTRIBUTION OF INTER-RANGE DISTANCES, THIS WOULD GUIDE MY CHOICE
# OF min.gapwidth.
# mtif.gaps = gaps(igt.mtifs.gr.rd)
# mtif.gaps.wth = width(mtif.gaps)
# hist(mtif.gaps.wth)
# 
# TOO BROAD A RANGE TO BE INFORMATIVE ABOUT THE SMALLEST GAPS
# mtif.gaps.wth.l100 = mtif.gaps.wth[which(mtif.gaps.wth <= 100)]
# library("ggplot2")
# mtif.gaps.wth.l100.df = as.data.frame(mtif.gaps.wth.l100)
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
# 
# PLOT THE SAME FOR THE WHOLE VECTOR FOR COMPARISON
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
# # [1]  1  5  5  7  9  9 10 16 18 19
# 
# EFFECT OF VARYING MIN.GAPWIDTH
# length(reduce(igt.mtifs.gr, drop.empty.ranges = FALSE, min.gapwidth = 0, with.revmap = TRUE))
# # [1] 2998
# length(reduce(igt.mtifs.gr, drop.empty.ranges = FALSE, min.gapwidth = 1, with.revmap = TRUE))
# # [1] 2993
# 
# DO THE SET OPERATION ON THE DATAFRAMES USING DPLYR - WHAT IS DIFFERENT
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
# SO IT SEEMS THAT THERE ARE 5 PAIRS OF ASSEMBLIES THAT ARE SUBSUMED BY min.gapwidth = 1
# THEREFORE, I THINK IT'S THE BEST THING TO TO STIPULATE min.gapwidth = 0
# 
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
# 
# WIDTH AND NO.READS DISTRIBUTION OF IGT MTIF ASSEMBLIES
library(ggplot2)
ggplot(igt.mtif.asm, aes(x = width)) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(0, 6000, by = 1000)) +  
  ylab("Counts") +
  xlab("IGT assembly width") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 26),
        axis.title.y = element_text(vjust = 1, size = 26),
        axis.text.x = element_text(size=24, vjust = 0.5),
        axis.text.y = element_text(size=24, vjust = 0.5),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
ggsave("igt-mtif-asm-wdth.png", dpi = 400)

ggplot(igt.mtif.asm, aes(x = log10(ypd.counts))) +
  geom_histogram() +
  ylab("Counts") +
  xlab("log10(reads in YPD)") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 26),
        axis.title.y = element_text(vjust = 1, size = 26),
        axis.text.x = element_text(size=24, vjust = 0.5),
        axis.text.y = element_text(size=24, vjust = 0.5),
        plot.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
ggsave("igt-mtif-asm-ypdReads.png", dpi = 400)
# 
# WRTIE OUT FOR POSTERITY
write.table(igt.mtif.asm,
          file = "igt-mtifs-asm.txt",
          sep = "\t",
          row.names = FALSE,
          quote = FALSE)
# 
# COULD DO ANOTHER VERSION OF THE GTF WITH NICK'S GTF WITH IGTS ATTACHED
# 
# REFORMATING INTO GTF FORMAT AND FUSING TO NICK'S GTF
system("awk -f anno-to-gtf2.awk igt-mtifs-asm-noHeader.txt > igt_mtif_asm.gtf")
system("cat sac_cer_yassour_utr.gtf igt_mtif_asm.gtf > sc_yassour_utr_steinmetz_igt.gtf")
# 
# NOW TO THE MORE USUAL MTIFS
library(plyr)
library(dplyr)
# 
coio.mtif.summ = coio.mtifs %>%
  group_by(gene.id) %>%
  summarise(sum.ypd.counts = sum(ypd.counts),
            no.mTifs = length(ypd.counts))
# 
# DISTRIBUTION OF NUMBER OF FORMS
ggplot(coio.mtif.summ, aes(x = log2(no.mTifs))) +
  geom_histogram() +
  ylab("Counts") +
  xlab("No. mTIFs covering one intact ORF per ORF") +
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
ggsave("no-mtifs-each-coio.png", dpi = 400)
# 
# FILTER EACH GENE, REMOVE THE BOTTOM 20% BY cume_dist()
coio.mtif.filt = coio.mtifs %>%
  group_by(gene.id) %>%
  filter(cume_dist(ypd.counts) > 0.2)
# 
length(unique(coio.mtif$gene.id)) == length(unique(coio.mtif.filt$gene.id))
# # [1] TRUE
length(unique(coio.mtifs$gene.id))
# # [1] 4729
dim(coio.mtifs)
# # [1] 136687      7
length(unique(coio.mtif.filt$gene.id))
# # [1] 4729
dim(coio.mtif.filt)
# [1] 123719      7
# 
# CONVERT FILTERED MTIFS INTO A GRANGES OBJECT
coio.mtif.filt.gr = makeGRangesFromDataFrame(coio.mtif.filt,
                                        keep.extra.columns = TRUE,
                                        ignore.strand = FALSE,
                                        seqnames.field = c("chr"),
                                        start.field = c("start"),
                                        end.field = c("end"),
                                        strand.field = c("strand"))
# 
# SPLIT INTO A GRANGESLIST BY GENE.ID
coio.mtif.filt.gr.l = split(coio.mtif.filt.gr, 
                            coio.mtif.filt.gr$gene.id, 
                            drop = TRUE)
# 
# SUBSUME THE ALIGNMENTS FOR EACH GENE BY LAPPLY-ING REDUCE()
coio.mtif.filt.gr.l.rd = lapply(coio.mtif.filt.gr.l, 
                                reduce, 
                                drop.empty.ranges = FALSE,
                                min.gapwidth = 0,
                                with.revmap = FALSE)
# 
# CAN'T SEEM TO UNLIST THIS, SO I'M GOING TO MAKE IT INTO A DATAFRAME, THEN CONVERT IS BACK INTO A GRANGES OBJECT
coio.mtif.filt.gr.l.rd = lapply(coio.mtif.filt.gr.l.rd, 
                                as.data.frame, 
                                row.names = NULL)
coio.mtif.filt.df = ldply(coio.mtif.filt.gr.l.rd)
colnames(coio.mtif.filt.df) = c("gene.id", "chr", "start", "end", "length", "strand")
# 
# WRITE OUT AS TABLE TO REFORMAT IN AWK
write.table(coio.mtif.filt.df,
            file = "coio-mtif-filt.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
# 
# REFORMAT WITH coio-mtif-filt-reformatter.awk
system('awk -f coio-mtif-filt-reformatter.awk coio-mtif-filt.txt > coio-mtif-filt.gtf')
# 
# JOIN WITH ONLY THE EXON LINES OF THE YASSOUR GTF FOR THE SUBSUMATION
# TAKE JUST THE LINES WITH EXON
system('grep -w "exon" sac_cer_yassour_utr.gtf | awk -f make-normal-gene-names.awk > sac_cer_yassour_exon.gtf')
# 
# JOIN USING CAT
system('cat coio-mtif-filt.gtf sac_cer_yassour_exon.gtf > sac_cer_exon.gtf')
#
# READ IN THE NEW GTF, SPLIT IT INTO A LIST BY GENE NAME, THE USED REDUCE TO SUBSUME OVERLAPPING
# ALIGNMENTS
library(rtracklayer)
exon.gtf = import("sac_cer_exon.gtf", format = "GFF", asRangedData = FALSE)
exon.gtf
length(exon.gtf)
names(exon.gtf) = mcols(exon.gtf)$group
exon.gtf.l = split(exon.gtf, exon.gtf$group)
length(unique(exon.gtf$group))
length(exon.gtf.l)
el.before = elementLengths(exon.gtf.l)
table(el.before)
# el.before
#    1    2    3    4    5    6    8 
# 1885 4640  157   11    1    1    1 
test4 = exon.gtf.l[which(elementLengths(exon.gtf.l) == 4)][[1]]
exon.gtf.l.rd = lapply(exon.gtf.l, 
                       reduce, 
                       drop.empty.ranges = FALSE,
                       min.gapwidth = 0,
                       with.revmap = FALSE)
el.after = elementLengths(exon.gtf.l.rd)
table(el.after)

elementLengths(exon.gtf.l[1])

foo = exon.gtf.l$YAL002W
sum(elementLengths(foo))

if (sum(elementLengths(foo)) == 2) {
  if (length(unique(foo$source))) {
    bar = subset(foo, foo$source == "steinmetz_mTIFs_coio")
  }
}
bar
rm(bar)
exon.gtf.l[[1]]$source
test3$source == "steinmetz_mTIFs_coio"
ranges(test3)
start(test3)[which(test3$source == "steinmetz_mTIFs_coio")]
end(test3)
max(start(test3)[which(test3$source == "utr-analysis")])

test4

for (i in 1:length(exon.gtf.l)) {
  if (sum(elementLengths(exon.gtf.l[[i]])) == 2) {
    if (length(unique(exon.gtf.l[[i]]$source)) == 2) {
      exon.gtf.l[[i]] = subset(exon.gtf.l[[i]], exon.gtf.l[[i]]$source == "steinmetz_mTIFs_coio")
    }
  } else if (sum(elementLengths(exon.gtf.l[[i]])) == 3) {
    if (length(unique(exon.gtf.l[[i]]$source)) == 2) {
      exon.gtf.l[[i]] = GRanges(seqnames = seqnames(exon.gtf.l[[i]])[1:length(exon.gtf.l[[i]])-1],
          ranges = IRanges(start = c(start(exon.gtf.l[[i]])[which(exon.gtf.l[[i]]$source == "steinmetz_mTIFs_coio")],
                                     max(start(exon.gtf.l[[i]])[which(exon.gtf.l[[i]]$source == "utr-analysis")])),
                           end = c(min(end(exon.gtf.l[[i]])[which(exon.gtf.l[[i]]$source == "utr-analysis")]),
                                   end(exon.gtf.l[[i]])[which(exon.gtf.l[[i]]$source == "steinmetz_mTIFs_coio")])),
          strand = strand(exon.gtf.l[[i]])[1:length(exon.gtf.l[[i]])-1],
          mcols = mcols(exon.gtf.l[[i]])[which(exon.gtf.l[[i]]$source == "utr-analysis"),])
      names(exon.gtf.l[[i]]) = mcols(exon.gtf.l[[i]])[,5]
      mcols(exon.gtf.l[[i]])$mcols.source = "steinmetz_mTIFS_coio"
    }
  }
}



coio.mtif.filt.rd.gr = makeGRangesFromDataFrame(coio.mtif.filt.df,
                                                keep.extra.columns = TRUE,
                                                ignore.strand = FALSE,
                                                seqnames.field = c("chr"),
                                                start.field = c("start"),
                                                end.field = c("end"),
                                                strand.field = c("strand"))
# 
# TAKE JUST THE LINES WITH EXON
system('grep -w "exon" sac_cer_yassour_utr.gtf > sac_cer_yassour_exon.gtf')
# 
# READ IN NICK'S GTF
library(rtracklayer)
gtf = import("sac_cer_yassour_exon.gtf", format = "GFF", asRangedData = FALSE)
gtf
coio.mtif.filt.rd.gr
# chr names are different
# gene names are not so accessible


