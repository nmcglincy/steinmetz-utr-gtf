# Attempt 3, using just the mTIFs covering one intact ORF, 20150211
# 
# Taking rows pertaining to TIFs covering one intact ORF, remove those with a gene name of 'NA',
# then reformat.
system("grep -w 'Covering one intact ORF' S2_tcd_mTIFAnno.txt | grep -wv 'NA' | awk -f reformater.awk > coio-mtifs.txt")
system("awk -f reorder-loc.awk coio-mtifs.txt > coio-mtifs2.txt")
# 
mtifs = read.delim("coio-mtifs2.txt",
                   header = FALSE)
colnames(mtifs) = c("chr", "strand", "start", "end", "ypd.counts", "class", "gene.id")
head(mtifs)
tail(mtifs)
# 
# FILTER EACH GENE, REMOVE THE BOTTOM 20% BY cume_dist()
library(plyr); library(dplyr)
mtifs.filt = mtifs %>%
  group_by(gene.id) %>%
  filter(cume_dist(ypd.counts) > 0.2)
# 
dim(mtifs)
dim(mtifs.filt)
length(unique(mtifs$gene.id))
length(unique(mtifs.filt$gene.id))
# 
# CONVERT FILTERED MTIFS INTO A GRANGES OBJECT
library("GenomicRanges")
mtifs.filt.gr = makeGRangesFromDataFrame(mtifs.filt,
                                         keep.extra.columns = TRUE,
                                         ignore.strand = FALSE, 
                                         seqnames.field = c("chr"),
                                         start.field = c("start"),
                                         end.field = c("end"),
                                         strand.field = c("strand"))
# 
# SPLIT INTO A GRANGESLIST BY GENE.ID
mtifs.filt.gr.l = split(mtifs.filt.gr, 
                        mtifs.filt.gr$gene.id, 
                        drop = TRUE)
# mtifs.filt.gr.l
# 
# SUBSUME THE ALIGNMENTS FOR EACH GENE BY LAPPLY-ING REDUCE()
mtifs.filt.gr.l.rd = lapply(mtifs.filt.gr.l, 
                            reduce, 
                            drop.empty.ranges = FALSE,
                            min.gapwidth = 0,
                            with.revmap = FALSE)
# mtifs.filt.gr.l.rd
# 
# CAN'T SEEM TO UNLIST THIS, SO I'M GOING TO MAKE IT INTO A DATAFRAME, THEN CONVERT IS BACK INTO A GRANGES OBJECT
mtifs.filt.gr.l.rd = lapply(mtifs.filt.gr.l.rd, 
                            as.data.frame, 
                            row.names = NULL)
mtif.filt.rd.df = ldply(mtifs.filt.gr.l.rd)
colnames(mtif.filt.rd.df) = c("gene.id", "chr", "start", "end", "length", "strand")
# 
# WRITE OUT AS TABLE TO REFORMAT IN AWK
write.table(mtif.filt.rd.df,
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
names(exon.gtf) = mcols(exon.gtf)$group
# 
head(exon.gtf)
tail(exon.gtf)
# 
exon.gtf.l = split(exon.gtf, exon.gtf$group)
# 
length(exon.gtf.l)
table(elementLengths(exon.gtf.l))
# 
source("EndPicker.R")
exon.gtf.l2 = lapply(exon.gtf.l, EndPicker)
table(elementLengths(exon.gtf.l2))

exon.gtf.l2[which(elementLengths(exon.gtf.l) > 5)]
exon.gtf.l2[which(elementLengths(exon.gtf.l2) > 5)]

length(exon.gtf.l2)
head(exon.gtf.l2)
names(exon.gtf.l2)

table(elementLengths(exon.gtf.l))
table(elementLengths(exon.gtf.l2))

exon.gtf.l[["YIL082W-A"]]
exon.gtf.l2[["YPR170W-B"]]

exon.gtf.gr = unlist(exon.gtf.l2,
                     recursive = TRUE,
                     use.names = FALSE)
class(exon.gtf.gr)
foo = unlist(exon.gtf.gr)
head(foo)
bar = lapply(foo, as.data.frame, row.names = NULL)
head(bar)
bar = lapply(bar, 
             function(x) {names(x) <- c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "phase", "group");
                          return(x)})
ape = ldply(bar)
head(ape)
length(unique(ape$group))
# this is the right length, 6696
#
# I want to check that there are no overlaps between genes on the same strand
# If this is true for a strand of x, then the length of disjoin(x) should be the same as the 
# length of x
# 
# convert ape df into GRanges
head(ape[,-4])
ape.gr = makeGRangesFromDataFrame(ape[,-4],
                                  keep.extra.columns = TRUE,
                                  ignore.strand = FALSE, 
                                  seqnames.field = c("seqnames"),
                                  start.field = c("start"),
                                  end.field = c("end"),
                                  strand.field = c("strand"))
head(ape.gr)
length(unique(ape.gr$group))
length(unique(ape$group))
dim(ape)
# 
?disjoin
monkey = disjoin(ape.gr, ignore.strand = FALSE)
head(monkey)
length(monkey)
# [1] 7790, well that's bad news
# banana = findOverlaps(ape.gr, ape.gr) obviously everything overlaps - not very useful
# 
library(rtracklayer)
export(ape.gr,
       con = "hybrid.gtf",
       format = "gtf")









