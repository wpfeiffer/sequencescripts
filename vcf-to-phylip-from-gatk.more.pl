#!/usr/bin/perl
use Getopt::Long;
use List::Util qw[min max];
$samples = 2;
$mixed = 0;
$more = 0;
$min_minor_pct = 5;
$min_minor_cov = 5;
GetOptions ("samples=i" => \$samples,
            "mixed=i" => \$mixed,
            "mixed_ids=s" => \$mixed_ids,
            "more=i" => \$more,
            "more_ids=s" => \$more_ids,
            "sep_more_pcts=s" => \$sep_more_pcts,
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
      if ($more > 0) {
        @more_id = split /,/, $more_ids;
        @sep_more_pct = split /,/, $sep_more_pcts;
      }
      for ($imix = 0; $imix < $mixed; $imix++) {
        $mixed_sites[$imix] = 0;
        $more_sites[$imix] = 0;
# Make room for mixed sites
        if ($mixed_id[$imix] ne $more_id[$imix]) {
          for ($j = $i; $sample_id[$j] ne $mixed_id[$imix]; $j--) {
            $sample_id[$j] = $sample_id[$j-1];
          }
          $sample_id[$j] = join ".", $sample_id[$j], "1";
          $sample_id[$j+1] = join ".", $sample_id[$j+1], "2";
# Make room for more sites
        } else {
          $i++;
          for ($j = $i; $sample_id[$j] ne $mixed_id[$imix]; $j--) {
            $sample_id[$j] = $sample_id[$j-2];
          }
          $sample_id[$j+1] = $sample_id[$j];
          $sample_id[$j] = join ".", $sample_id[$j], "1";
          $sample_id[$j+1] = join ".", $sample_id[$j+1], "2";
          $sample_id[$j+2] = join ".", $sample_id[$j+2], "3";
        }
        $i++;
      }
    }
  }
# Consider only SNPs
  if (!($fields[0] =~ /#/) && length($fields[3]) == 1 && length($fields[4]) == 1) {
    $ref = $fields[3];
    $var = $fields[4];
    for ($i = 0; $i < $samples; $i++) {
      $sample_field = 9 + $i;
      ($gt, $cov) = split /:/, $fields[$sample_field];
      ($ref_cov[$i], $var_cov[$i]) = split /,/, $cov;
# Overrule genotype from GATK
      if ($ref_cov[$i] >= $var_cov[$i]) {
        $phylip[$i][$nvar] = $ref;
      } else {
        $phylip[$i][$nvar] = $var;
      }
    }
    if ($mixed > 0) {
      for ($imix = 0; $imix < $mixed; $imix++) {
        $mixed_idjoin = join ".", $mixed_id[$imix], "1";
# Initialize mixed sites
        if ($mixed_id[$imix] ne $more_id[$imix]) {
          for ($j = $i; $sample_id[$j] ne $mixed_idjoin; $j--) {
            $phylip[$j][$nvar] = $phylip[$j-1][$nvar];
            $ref_cov[$j] = $ref_cov[$j-1];
            $var_cov[$j] = $var_cov[$j-1];
          }
# Initialize more sites
        } else {
          $i++;
          for ($j = $i; $sample_id[$j] ne $mixed_idjoin; $j--) {
            $phylip[$j][$nvar] = $phylip[$j-2][$nvar];
            $ref_cov[$j] = $ref_cov[$j-2];
            $var_cov[$j] = $var_cov[$j-2];
          }
          $phylip[$j+1][$nvar] = $phylip[$j][$nvar];
        }
        $i++;
        $tot_cov = max(1, $ref_cov[$j]+$var_cov[$j]);
        if (100*$ref_cov[$j]/$tot_cov > $min_minor_pct && 100*$var_cov[$j]/$tot_cov > $min_minor_pct && $ref_cov[$j] >= $min_minor_cov && $var_cov[$j] >= $min_minor_cov) {
# Evaluate mixed variant sites
          if ($mixed_id[$imix] ne $more_id[$imix]) {
            $mixed_sites[$imix]++;
            if ($ref_cov[$j] >= $var_cov[$j]) {
              $phylip[$j][$nvar] = $ref;
              $phylip[$j+1][$nvar] = $var;
              $minor_gt = "var";
              $minor_pct = 100*$var_cov[$j]/$tot_cov;
# Evaluate mixed reference sites
            } else {
              $phylip[$j][$nvar] = $var;
              $phylip[$j+1][$nvar] = $ref;
              $minor_gt = "ref";
              $minor_pct = 100*$ref_cov[$j]/$tot_cov;
            }
            ($dum1, $dum2, $dum3, $contig_id) = split /\Q|/, $fields[0];
            printf MIXFILE "%12s  %14s  %8d  %3s  %5.2f\n", $sample_id[$j+1], $contig_id, $fields[1], $minor_gt, $minor_pct;
# Evaluate more variant sites
          } else {
            if ($ref_cov[$j] >= $var_cov [$j]) {
              $phylip[$j][$nvar] = $ref;
# If $sep_more_pct > 0, then more sites are independent of mixed sites
              if ($sep_more_pct[$imix] > 0) {
                if (100*$var_cov[$j]/$tot_cov > abs($sep_more_pct[$imix])) {
                  $mixed_sites[$imix]++;
                  $phylip[$j+1][$nvar] = $var;
                  $minor_gt = "var";
                  $minor_pct = 100*$var_cov[$j]/$tot_cov;
                  $phylip[$j+2][$nvar] = $ref;
                  $more_pct = 0;
                } else {
                  $phylip[$j+1][$nvar] = $ref;
                  $minor_pct = 0;
                  $more_sites[$imix]++;
                  $phylip[$j+2][$nvar] = $var;
                  $more_gt = "var";
                  $more_pct = 100*$var_cov[$j]/$tot_cov;
                }
# If $sep_more_pct < 0, then more sites are subsidiary to mixed sites
              } else {
                if (100*$var_cov[$j]/$tot_cov > abs($sep_more_pct[$imix])) {
                  $mixed_sites[$imix]++;
                  $phylip[$j+1][$nvar] = $var;
                  $minor_gt = "var";
                  $minor_pct = 100*$var_cov[$j]/$tot_cov;
                  $more_sites[$imix]++;
                  $phylip[$j+2][$nvar] = $var;
                  $more_gt = "var";
                  $more_pct = 100*$var_cov[$j]/$tot_cov;
                } else {
                  $phylip[$j+1][$nvar] = $ref;
                  $minor_pct = 0;
                  $more_sites[$imix]++;
                  $phylip[$j+2][$nvar] = $var;
                  $more_gt = "var";
                  $more_pct = 100*$var_cov[$j]/$tot_cov;
                }
              }
# Evaluate more reference sites
            } else {
              $phylip[$j][$nvar] = $var;
# If $sep_more_pct > 0, then more sites are independent of mixed sites
              if ($sep_more_pct[$imix] > 0) {
                if (100*$ref_cov[$j]/$tot_cov > abs($sep_more_pct[$imix])) {
                  $mixed_sites[$imix]++;
                  $phylip[$j+1][$nvar] = $ref;
                  $minor_gt = "ref";
                  $minor_pct = 100*$ref_cov[$j]/$tot_cov;
                  $phylip[$j+2][$nvar] = $var;
                  $more_pct = 0;
                } else {
                  $phylip[$j+1][$nvar] = $var;
                  $minor_pct = 0;
                  $more_sites[$imix]++;
                  $phylip[$j+2][$nvar] = $ref;
                  $more_gt = "ref";
                  $more_pct = 100*$ref_cov[$j]/$tot_cov;
                }
# If $sep_more_pct < 0, then more sites are subsidiary to mixed sites
              } else {
                if (100*$ref_cov[$j]/$tot_cov > abs($sep_more_pct[$imix])) {
                  $mixed_sites[$imix]++;
                  $phylip[$j+1][$nvar] = $ref;
                  $minor_gt = "ref";
                  $minor_pct = 100*$ref_cov[$j]/$tot_cov;
                  $more_sites[$imix]++;
                  $phylip[$j+2][$nvar] = $ref;
                  $more_gt = "ref";
                  $more_pct = 100*$ref_cov[$j]/$tot_cov;
                } else {
                  $phylip[$j+1][$nvar] = $var;
                  $minor_pct = 0;
                  $more_sites[$imix]++;
                  $phylip[$j+2][$nvar] = $ref;
                  $more_gt = "ref";
                  $more_pct = 100*$ref_cov[$j]/$tot_cov;
                }
              }
            }
            ($dum1, $dum2, $dum3, $contig_id) = split /\Q|/, $fields[0];
            if ($minor_pct > 0) {
              printf MIXFILE "%12s  %14s  %8d  %3s  %5.2f\n", $sample_id[$j+1], $contig_id, $fields[1], $minor_gt, $minor_pct;
            } 
            if ($more_pct > 0) {
              printf MIXFILE "%12s  %14s  %8d  %3s  %5.2f\n", $sample_id[$j+2], $contig_id, $fields[1], $more_gt, $more_pct;
            }
          }
        }
      }
    }
    $nvar++;
  }
}
if ($nvar > 0) {
  $taxa = $samples + $mixed + $more;
  print "$taxa $nvar\n";
  for ($i = 0; $i < $taxa; $i++) {
    print "$sample_id[$i] ", @{$phylip[$i]}, "\n";
  }
  if ($mixed > 0) {
    print MIXFILE "min_minor_pct $min_minor_pct\n";
    print MIXFILE "min_minor_cov $min_minor_cov\n";
    print MIXFILE "mixed $mixed\n";
    print MIXFILE "mixed_ids @mixed_id\n";
    print MIXFILE "mixed_sites @mixed_sites\n";
  }
  if ($more > 0) {
    print MIXFILE "more $more\n";
    print MIXFILE "more_ids @more_id\n";
    print MIXFILE "sep_more_pcts @sep_more_pct\n";
    print MIXFILE "more_sites @more_sites\n";
  }
}
