#!/bin/awk -f
# AN AWK SCRIPT TO CHANGE LOCATION DETAILS TO MAKE IT CONSISTENT WITH NICK'S GTF.
# SMALLER COORDINATE NEEDS TO BE FIRST EVEN IF IT'S ON THE NEGATIVE STRAND
# 20150102 - ALSO REMOVES THE 6TH FIELD, WHICH SHOULD BE THE YGAL COUNTS
BEGIN {OFS = "\t"}
{
    if ( $5 >= 1 ) {
        if ( $2 == "-") {
            print $1, $2, $4, $3, $5, $7, $8
        } else {
            print $1, $2, $3, $4, $5, $7, $8
        }
    }
}