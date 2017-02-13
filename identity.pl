#!/usr/bin/perl
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
