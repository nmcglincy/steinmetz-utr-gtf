#!/bin/awk -f
# A SCRIPT TO CONVERT MY POLISHED ANNOTATION TO A REAL GTF
# A GTF SHOULD HAVE THE FOLLOWING FIELDS
# 1 <seqname> - generally the chromosome
# 2 <source> 
# 3 <feature> - should all be 'exon'
# 4 <start>
# 5 <end>
# 6 <score> - use number of reads here
# 7 <strand>
# 8 <frame> - nick puts '.' here
# 9 [attributes] [comments] - ; separated list
BEGIN { OFS = "\t"}
{
    if ( length($1) < 2) {
        print "chr0" $1 OFS "steinmetz_mTIFs" OFS "exon" OFS $3 OFS $4 OFS $5 OFS $2 OFS "." OFS "gene_id  \"" $7 "\"; transcript_id \"" $7 "\"; mTIF_type \"" $6 "\";"
    } else {
        print "chr" $1 OFS "steinmetz_mTIFs" OFS "exon" OFS $3 OFS $4 OFS $5 OFS $2 OFS "." OFS "gene_id  \"" $7 "\"; transcript_id \"" $7 "\"; mTIF_type \"" $6 "\";"
    }
}