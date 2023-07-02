use strict;
use warnings;
use List::Util qw[min max];
my @hit;
my $i;
my $j;
my $k;

my $usage = "Usage: cat <outfile> | perl framemark-coverage.pl <posfile> \n";

my $posfile = shift;
my $eval_cuttoff = shift;

my @fields;
my %pos_starts = ();
my %pos_ends = ();
my $fam;
my $start;
my $end;
my $line;
my @temp;

#get all positive hits from outfiles
$i = 0;
while($line = <>) { 
    chomp $line;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    @temp = split(/\s+/, $line);
    # 0      1          2             3           4        5       6              7
    #-8     -7         -6            -5          -4       -3      -2             -1
    #($pval, $bitscore, $target_from, $target_to, $target, $query, $matching_hit, $matching_strand) 
    # 0      1          2             3           4        5       6    7      8      9               10
    #-11    -10        -9            -8          -7       -6      -5   -4     -3     -2              -1
    #($pval, $bitscore, $target_from, $target_to, $target, $query, $fs, $stop, $pipe, $matching_hit,  $matching_strand)
    unless($temp[-2]  =~ m/^decoy/ || ($temp[-2] !~ m/^$temp[5]\/\d+/) || ($temp[-1] !~ /same/)) {   
    	push(@{$hit[$i]}, split(/\s+/, $line)); 
    	$i++;
    }
}

#get true coordinates of each embemdded sequence 
open(POS, $posfile) || die "ERROR, couldn't open $posfile";
while ($line = <POS>)
{
    chomp $line;
    if ($line =~ m/^\#/) { next; }
    @fields   = split(' ', $line, 5);
    $fam = $fields[0];
    $start = $fields[3];
    $end   = $fields[4];
    $pos_starts{$fam} = min($start, $end);
    $pos_ends{$fam}    = max($start, $end);
}

my %matches = ();
my $name;
my @froms;
my @tos;
my $from;
my $to;
my $s;
my $e;
my $eval_over;
my $long_over;
my $total_over;
my $eval_cover;
my $long_cover;
my $total_cover;
my $tmp_long;
my $tmp_over;
my $true_length;

printf("%20s	%11s	%10s	%10s	%10s	%10s	%10s	%10s	%10s\n", "HIT",  "TRUE-LENGTH", "EVAl-COVER", "EVAL-OVER", "LONG-COVER", "LONG-OVER", "TOT-COVER", "TOT-OVER", "CONTIGS");
for($i = 0; $i < @hit; $i++) {
	$name = $hit[$i][-2];	
    	unless(exists($matches{$name})) {
		if($hit[$i][0] <= $eval_cuttoff) {
			$matches{$name} = 1;
			@froms          = ();
			@tos            = ();
			$from           = $hit[$i][2];
                        $to             = $hit[$i][3];
			push(@froms, min($from,$to));
			push(@tos,   max($from,$to));

			#hits are already sorted by evalue so first hit has lowest evalue
			$eval_cover  = $tos[0] - $froms[0] + 1;
			$eval_over   = max($pos_starts{$name} - $froms[0], 0);
			$eval_over  += max($pos_ends{$name}   - $tos[0],   0);
			$true_length  = $pos_ends{$name} - $pos_starts{$name} + 1; 

			for($j = $i+1; $j < @hit; $j++) {
				if($hit[$j][6] eq $name && $hit[$j][0] <= 0.00001) {
					$matches{$name} += 1;
					$from = $hit[$j][2];
					$to   = $hit[$j][3];
					push(@froms, min($from,$to));
					push(@tos,   max($from,$to));				
				}
			}
	
			#if only on hit it is the longest
			if($matches{$name} == 1) {
                                $long_cover  = $tos[0] - $froms[0] + 1;
				$long_over   = $pos_starts{$name} - $froms[0];
				$long_over  += $pos_ends{$name}   - $tos[0];
				$long_cover -= $long_over;
                        }
			#if multiple hits find the hit with best coverage
			else {
				$long_cover  = $tos[0] - $froms[0] + 1;
				$long_over   = $pos_starts{$name} - $froms[0];
				$long_over  += $pos_ends{$name}   - $tos[0];
				$long_cover -= $long_over;

				for($k = 1; $k <  @froms; $k++) {
					if(($tos[$k] - $froms[$k] + 1) > $long_cover) { #new hit is longer so might have better coverage
						$tmp_long  = $tos[$k] - $froms[$k] + 1;
						$tmp_over  = $pos_starts{$name} - $froms[$k];
						$tmp_over += $pos_ends{$name}   - $tos[$k];
						$tmp_long -= $tmp_over;
						if($tmp_long > $long_cover) { #new hit has better coverage
							$long_cover = $tmp_long;
							$long_over  = $tmp_over;
						} 
						elsif($tmp_long == $long_cover) { #new hit has equal coverage pick the one with lower overage
							if($long_over > $tmp_over) { #new hit has lover overage
								$long_cover = $tmp_long;
								$long_over  = $tmp_over;
							}
						}					
					}
				}
			}

			#find total overage														
			$s          = min(@froms);
                        $e          = max(@tos);
			$total_over = 0;

                        if($s < $pos_starts{$name}) {
                                $total_over += $pos_starts{$name} - $s;
                                $s           = $pos_starts{$name};
                        }
                        if($e > $pos_ends{$name}) {
                                $total_over += $e - $pos_ends{$name};
                                $e           = $pos_ends{$name};
                        }

			#sort hits by start 
			my @idx = sort { $froms[$a] <=> $froms[$b] } 0 .. $#froms;
			@froms = @froms[@idx];
			@tos = @tos[@idx];    	
		
			#find total coverage
			$total_cover = 0;

			#if only one hit use its coverage as total
			if($matches{$name} == 1) {
				$total_cover = $e - $s + 1;
			}
			#if multiple hits consider all hits	
			else {
				#first remove hits which are entrily overlapped by other hits
				#they cannot add to coverage and just complcate calcualtions
				for($k = 1; $k <  @froms; $k++) {
					if($froms[$k] >= $froms[$k-1] and $tos[$k] <= $tos[$k+1]) {
						splice @froms, $k, 1;
						splice @tos,   $k, 1;
						$matches{$name} -= 1;
					}
					elsif($froms[$k-1] >= $froms[$k] and $tos[$k-1] <= $tos[$k]) {
						splice @froms, $k-1, 1;
						splice @tos,   $k-1, 1;
						$matches{$name}   -= 1;	
					}
				}
				#start with true hit length at total coverage then subtract all gaps
				$total_cover  = $pos_ends{$name} - $pos_starts{$name} + 1;

				#subtract gap any gap at begining and end using min start and 
				#max end variables from total overage calcuation
				$total_cover -= $s - $pos_starts{$name}; 
				$total_cover -= $pos_ends{$name} - $e;
			
				#if there are still multiple hits after removing total overlaps
				#look for geps between hits and subtract from total coverage
				for($k = 1; $k <  @froms; $k++) {
					if($tos[$k-1] < $froms[$k]) {
						$total_cover -= $froms[$k] - $tos[$k-1];
					}	
				}			
			}

			printf("%20s	%10d	%10d	%10d	%10d    %10d	%10d    %10d	%10d\n", $name, $true_length, $eval_cover, $eval_over, $long_cover, $long_over, $total_cover, $total_over, $matches{$name});
		}
	}
}

