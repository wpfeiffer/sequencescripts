#!/usr/bin/perl
$_ = <>;
print;
@fields = split;
$nsamp = $fields[0];
$nchar = $fields[1];
for ($i = 0; $i < $nsamp; $i++) {
  $_ = <>;
  @fields = split;
  $sample_id[$i] = $fields[0];
  @temp = split //, $fields[1];
  for ($k = 0; $k < $nchar; $k++) {
    $seq[$i][$k] = $temp[$k];
  }
}
for ($i = 0; $i < $nsamp; $i++) {
  for ($j = 0; $j < $nsamp; $j++) {
    $dist[$i][$j] = 0;
    if ($j > $i) {
      for ($k = 0; $k < $nchar; $k++) {
        if ($seq[$i][$k] ne $seq[$j][$k] && $seq[$i][$k] ne "-" && $seq[$j][$k] ne "-") { $dist[$i][$j]++; } 
      }
    }
    if ($j < $i) { $dist[$i][$j] = $dist[$j][$i]; }
  }
}
for ($i = 0; $i < $nsamp; $i++) {
  printf "%-20s" . ("%4d" x $nsamp) . "\n", $sample_id[$i], @{$dist[$i]};
}
