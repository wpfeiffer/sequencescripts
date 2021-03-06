#!/usr/bin/perl
# extract-high-confidence-from-gatk picks variants deemed of high confidence from the vcf file generated by GATK and outputs the variants to another vcf file.
# Such variants have each allele represented by at least one sample with
# (1) the percentage of reads corresponding to the minor allele is no more than max_minor_pct, and
# (2) the read coverage is between min_cov and max_cov inclusive.
# In addition,
# (3) the variant cannot be within min_space sites of another high-confidence variant for either of them to pass,
# (4) all samples need to have coverage ≥ min_cov if all_high_cov=yes, and
# (5) no more than nmix samples can have mixed variants, i.e., ones for which the percentage of reads corresponding to each allele is more than max_minor_pct.
#
# Input parameters are as follows:
# samples is the number of samples.
# min_space is the minimum spacing (+1) between a pair of variant sites for the variants to pass.
# min_cov is the minimum read coverage at the variant site for the variant to pass, unless overridden by all_high_cov=no.
# max_cov is the maximum read coverage at the variant site for the variant to pass.
# max_minor_pct is the maximum percentage of reads in the minor allele for a sample not to be considered mixed.
# max_mix is the maximume number of samples with mixed alleles at the variant site for the variant to pass.
# all_high_cov is a yes/no flag that specifies whether all samples need to have high coverage for the variant to pass.
#
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
