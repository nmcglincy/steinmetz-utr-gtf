#!/bin/awk -f
#
# AN AWK SCRIPT TO REFORMAT SAC_CER_YASSOUR.GTF GENE NAMES INTO SOMETHING EASIER TO DEAL WITH
#
BEGIN { FS = "\t|\""; OFS = "\t" }
{
	print $1, $2, $3, $4, $5, $6, $7, $8, $10
}