EndPicker = function(x){
  if (sum(elementLengths(x)) == 1) {
    x = x
  } else if (sum(elementLengths(x)) == 2) {
    if (length(unique(x$source)) == 2) {
      x = subset(x, x$source == "steinmetz_mTIFs_coio")
    } else {
      x = x
    }
  } else if (sum(elementLengths(x)) == 3) {
    if (length(unique(x$source)) == 2) {
      x = GRanges(seqnames = seqnames(x)[1:length(x)-1],
                                ranges = IRanges(start = c(start(x)[which(x$source == "steinmetz_mTIFs_coio")],
                                                         max(start(x)[which(x$source == "utr-analysis")])),
                                                 end   = c(min(end(x)[which(x$source == "utr-analysis")]),
                                                         end(x)[which(x$source == "steinmetz_mTIFs_coio")])),
                                strand = strand(x)[1:length(x)-1],
                                mcols  = mcols(x)[which(x$source == "utr-analysis"),])
    # names(x) = mcols(x)[,5]
    # mcols(x)$mcols.source = "steinmetz_mTIFS_coio"
  } else {
    x = x
  }
} else if (sum(elementLengths(x)) > 3) {
  if (length(unique(x$source)) == 2) {
    x = GRanges(seqnames = seqnames(x)[1:length(x)-1],
                              ranges   = IRanges(start = c(start(x)[which(x$source == "steinmetz_mTIFs_coio")],
                                                           sort(start(x)[which(x$source == "utr-analysis")])[-1]),
                                                 end   = c(sort(end(x)[which(x$source == "utr-analysis")])[1:length(sort(end(x)[which(x$source == "utr-analysis")]))-1],
                                                           end(x)[which(x$source == "steinmetz_mTIFs_coio")])),
                              strand   = strand(x)[1:length(x)-1],
                              mcols    = mcols(x)[which(x$source == "utr-analysis"),])
    # names(x) = mcols(x)[,5]
    # mcols(x)$mcols.source = "steinmetz_mTIFS_coio"  
    } else {
      x = x
    }
  }
}