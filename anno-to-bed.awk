#!/bin/awk -f
# A SCRIPT TO CONVERT MY POLISHED ANNOTATION TO A REAL BED FILE
BEGIN { OFS = "\t"}
{
    if ( length($1) < 2) {
        print "chr0" $1 OFS $3 OFS $4 OFS $7 OFS $5 OFS $2
    } else {
        print "chr" $1 OFS $3 OFS $4 OFS $7 OFS $5 OFS $2
    }
}

#    V1       V2       V3           V4 V5 V6
#1 chr1 66999824 67210768    NM_032291  0  +
#2 chr1 33546713 33585995    NM_052998  0  +
#3 chr1 48998526 50489626    NM_032785  0  -
#4 chr1 16767166 16786584 NM_001145278  0  +
#5 chr1 16767166 16786584 NM_001145277  0  +
#6 chr1  8384389  8404227 NM_001080397  0  +