#!/usr/bin/perl
use Getopt::Long;
use List::Util qw[min max];
$samples = 2;
$min_cov = 20;
$max_cov = 1000;
$mixed = 0;
$min_minor_pct = 5;
$min_minor_cov = 5;
GetOptions ("samples=i" => \$samples,
            "min_cov=i" => \$min_cov,
            "max_cov=i" => \$max_cov,
            "mixed=i" => \$mixed,
            "mixed_ids=s" => \$mixed_ids,
            "mixed_file=s" => \$mixed_file,
            "min_minor_pct=i" => \$min_minor_pct,
            "min_minor_cov=i" => \$min_minor_cov);
open MIXFILE, ">$mixed_file";
$nvar = 0;
while (<>) {
  @fields = split;
  if ($fields[0] =~ /#C/) {
    for ($i = 0; $i < $samples; $i++) {
      $sample_field = 9 + $i;
      $sample_id[$i] = $fields[$sample_field];
    }
    if ($mixed > 0) {
      @mixed_id = split /,/, $mixed_ids;
      for ($imix = 0; $imix < $mixed; $imix++) {
        $mixed_sites[$imix] = 0;
        for ($j = $i - 1; $j >= 0; $j--) {
          if ($sample_id[$j] eq $mixed_id[$imix]) {
            $imixed_id[$imix] = $j;
          }
        }
        for ($j = $i; $j > $imixed_id[$imix]; $j--) {
          $sample_id[$j] = $sample_id[$j-1];
        }
        $sample_id[$j] = join ".", $sample_id[$j], "1";
        $sample_id[$j+1] = join ".", $sample_id[$j+1], "2";
        $i++;
      }
    }
  }
# Consider only SNPs
  if (!($fields[0] =~ /#/) && length($fields[3]) == 1 && length($fields[4]) == 1 ) {
    $ref = $fields[3];
    $var = $fields[4];
    for ($i = 0; $i < $samples; $i++) {
      $sample_field = 9 + $i;
      ($gt, $cov) = split /:/, $fields[$sample_field];
      ($ref_cov, $var_cov) = split /,/, $cov;
      $tot_cov = max(1, $ref_cov+$var_cov);
# overrule genotype from GATK
      if ($ref_cov >= $var_cov) {
        $phylip[$i][$nvar] = $ref;
      } else {
        $phylip[$i][$nvar] = $var;
      }
      if ($tot_cov < $min_cov || $tot_cov > $max_cov) {
        $phylip[$i][$nvar] = "-";
      }
      if (100*$ref_cov/$tot_cov > $min_minor_pct && 100*$var_cov/$tot_cov > $min_minor_pct) {
        $phylip[$i][$nvar] = "-";
      }
    }
    if ($mixed > 0) {
      for ($imix = 0; $imix < $mixed; $imix++) {
        for ($j = $i; $j > $imixed_id[$imix]; $j--) {
          $phylip[$j][$nvar] = $phylip[$j-1][$nvar];
        }
        $i++;
        $sample_field = 9 + $j - $imix;
        ($gt, $cov) = split /:/, $fields[$sample_field];
        ($ref_cov, $var_cov) = split /,/, $cov;
        $tot_cov = max(1, $ref_cov+$var_cov);
        if (100*$ref_cov/$tot_cov > $min_minor_pct && 100*$var_cov/$tot_cov > $min_minor_pct && $ref_cov >= $min_minor_cov && $var_cov >= $min_minor_cov && $tot_cov >= $min_cov && $tot_cov <= $max_cov) {
          $mixed_sites[$imix]++;
          if ($ref_cov >= $var_cov) {
            $phylip[$j][$nvar] = $ref;
            $phylip[$j+1][$nvar] = $var;
            $minor_gt = "var";
            $minor_pct = 100*$var_cov/$tot_cov;
          } else {
            $phylip[$j][$nvar] = $var;
            $phylip[$j+1][$nvar] = $ref;
            $minor_gt = "ref";
            $minor_pct = 100*$ref_cov/$tot_cov;
          }
          ($dum1, $dum2, $dum3, $contig_id) = split /\Q|/, $fields[0];
          printf MIXFILE "%12s  %14s  %8d  %3s  %5.2f\n", $sample_id[$j+1], $contig_id, $fields[1], $minor_gt, $minor_pct;
        }
      }
    }
#   $phylip[$i][$nvar] = $ref;
    $nvar++;
  }
}
if ($nvar > 0) {
  $taxa = $samples + $mixed;
  print "$taxa $nvar\n";
  for ($i = 0; $i < $taxa; $i++) {
    print "$sample_id[$i] ", @{$phylip[$i]}, "\n";
  }
# print "reference ", @{$phylip[$i]}, "\n";
  print MIXFILE "min_minor_pct $min_minor_pct\n";
  print MIXFILE "min_minor_cov $min_minor_cov\n";
  print MIXFILE "mixed $mixed\n";
  print MIXFILE "mixed_ids @mixed_id\n";
  print MIXFILE "mixed_sites @mixed_sites\n";
}
