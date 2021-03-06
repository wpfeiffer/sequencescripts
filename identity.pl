#!/usr/bin/perl
# identity calculates the average nucleotide identity (ANI) and total length of matches between a sample fasta file and a reference fasta file based upon a prior BLAST alignment.
#
# The BLAST output, which is input to identity, is assumed to have been generated by the following command:
# 
# blastn -query <sample>.fasta -subject <reference>.fasta -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send' -out <sample>.vs.<reference>.out
#
# Input parameters are as follows:
# min_length is the minimum length of matches to consider.
# min_pident is a threhold value in percent below which matches are neglected.
#
use Getopt::Long;
$min_length = 0;
$min_pident = 90;
GetOptions ("min_length=i" => \$min_length,
            "min_pident=i" => \$min_pident);
$avg_identity = 0;
$tot_length= 0;
while (<>) {
  @fields = split;
  $identity = $fields[2];
  $length = $fields[3];
  if ($length >= $min_length && $identity >= $min_pident) {
    $avg_identity += $identity*$length;
    $tot_length += $length;
  }
}
$avg_identity = $avg_identity/$tot_length;
printf " min_pident = %2d; avg_identity = %6.2f%%; tot_length = %7d bp\n", $min_pident, $avg_identity, $tot_length;
