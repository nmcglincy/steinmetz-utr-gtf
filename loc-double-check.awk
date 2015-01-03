#!/bin/awk -f
# AN AWK SCRIPT TO ENSURE THE GENOMIC COORDINATES ARE IN THE RIGHT ORDER
{
    if ($4 > $3) {
        print "looks good"
    } else {
        print "oh shit"
    }
    print
}