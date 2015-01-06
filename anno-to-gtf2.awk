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
BEGIN { OFS = "\t" }
{
    if ( length($1) < 2) {
        print "chr0" $1, "steinmetz_igt_mTIF_asm", "exon" , $2 , $3 , $7 , $5 , "." , "gene_id  \"" $6 "\"; transcript_id \"" $6 "\";"
    } else {
        print "chr" $1, "steinmetz_igt_mTIF_asm", "exon" , $2 , $3 , $7 , $5 , "." , "gene_id  \"" $6 "\"; transcript_id \"" $6 "\";"
    }
}
