#!/bin/awk -f
#
# AN AWK SCRIPT TO REFORMAT COIO-MTIF-FILT.TXT INTO A GTF FORMAT
#
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
#
BEGIN { OFS = "\t" }
{
	if ( length($2) < 2) {
        print "chr0" $2, "steinmetz_mTIFs_coio" , "exon" , $3 , $4 , "0" , $6 , "." , $1
    } else {
        print "chr" $2, "steinmetz_mTIFs_coio" , "exon" , $3 , $4 , "0" , $6 , "." , $1
    }
	
}
#
# OLD LAST FIELD
# "gene_id  \"" $1 "\"; transcript_id \"" $1 "\";"