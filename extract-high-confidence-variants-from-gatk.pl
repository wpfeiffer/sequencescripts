#!/usr/bin/perl
use Getopt::Long;
use List::Util qw[min max];
$samples = 2;
$min_space = 30;
$min_cov = 20;
$max_cov = 1000;
$max_minor_pct = 5;
$max_mix = 100;
$all_high_cov = "no";
GetOptions ("samples=i" => \$samples,
            "min_space=i" => \$min_space,
            "min_cov=i" => \$min_cov,
            "max_cov=i" => \$max_cov,
            "max_minor_pct=i" => \$max_minor_pct,
            "max_mix=i" => \$max_mix,
            "all_high_cov=s" => \$all_high_cov);
$l = 0;
while (<>) {
  @fields = split;
  if ($fields[0] =~ /#/) {
    print "$_";
  } else {
    $imin_cov = 1;
    $iprint_ref = 0;
    $iprint_var = 0;
    $nmix = 0;
    for ($i = 0; $i < $samples; $i++) {
      $sample_field = 9 + $i;
      ($gt[$i], $cov) = split /:/, $fields[$sample_field];
      ($ref_cov, $var_cov) = split /,/, $cov;
# overrule genotype from GATK
      if ($ref_cov >= $var_cov) {
        $gt[$i] = "0";
      } else {
        $gt[$i] = "1";
      }
      $tot_cov = max(1, $ref_cov+$var_cov);
# require all samples to have high coverage if $all_high_cov eq "yes"
      if ($tot_cov < $min_cov && $all_high_cov eq "yes") { $imin_cov = 0; }
# require high confidence and high coverage (but not too high) for at least one reference sample and one variant sample
      if ($gt[$i] eq "0" && 100*$var_cov/$tot_cov <= $max_minor_pct && $tot_cov >= $min_cov && $tot_cov <= $max_cov) { $iprint_ref = 1; }
      if ($gt[$i] eq "1" && 100*$ref_cov/$tot_cov <= $max_minor_pct && $tot_cov >= $min_cov && $tot_cov <= $max_cov) { $iprint_var = 1; }
# count sites with mixed alleles
      if (100*$ref_cov/$tot_cov > $max_minor_pct && 100*$var_cov/$tot_cov > $max_minor_pct) { $nmix++ }
    }
# check whether high confidence and high coverage variants are close
    if ($iprint_ref == 1 && $iprint_var == 1 && $imin_cov == 1 && $nmix <= $max_mix) {
      $line[$l] = $_;
      $chrom[$l] = $fields[0];
      $pos[$l] = $fields[1];
      $iclose[$l] = 0;
      if ($l > 0 && $chrom[$l] eq $chrom[$l-1] && ($pos[$l] - $pos[$l-1]) <= $min_space) {
        $iclose[$l-1] = 1;
        $iclose[$l] = 1;
      }
      $l++;
    }
  }
}
$nline = $l;
# print out only isolated variants
for ($l = 0; $l < $nline; $l++) {
  if ($iclose[$l] == 0) { print "$line[$l]"; }
}
