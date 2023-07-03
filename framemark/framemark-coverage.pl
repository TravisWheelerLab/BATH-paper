use strict;
use warnings;
use List::Util qw[min max];
my @hit;
my $i;
my $j;
my $k;

my $usage = "Usage: cat <outfile> | perl framemark-coverage.pl <posfile> <eval>\n";

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
my $total_over;
my $total_cover;
my $true_length;

printf("%20s	%11s	%10s	%10s	%10s\n", "# HIT",  "TRUE-LENGTH", "COVERED", "OVER", "CONTIGS");
for($i = 0; $i < @hit; $i++) {
	$name = $hit[$i][-2];	
	#if this is a new hit
    	unless(exists($matches{$name})) {
		#if the hit passes the eval cuttoff
		if($hit[$i][0] <= $eval_cuttoff) {
			$matches{$name} = 1;
			@froms          = ();
			@tos            = ();
			$from           = $hit[$i][2];
                        $to             = $hit[$i][3];
			push(@froms, min($from,$to));
			push(@tos,   max($from,$to));

			$true_length  = $pos_ends{$name} - $pos_starts{$name} + 1; 

			for($j = $i+1; $j < @hit; $j++) {
				if($hit[$j][6] eq $name && $hit[$j][0] <= $eval_cuttoff) {
					$matches{$name} += 1;
					$from = $hit[$j][2];
					$to   = $hit[$j][3];
					push(@froms, min($from,$to));
					push(@tos,   max($from,$to));				
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
				#first remove hits which are entrily overlapped by a previous hit
				#they cannot add to coverage and just complcate calcualtions
				for($k = 1; $k <  @froms; $k++) {
					if($tos[$k] <= $tos[$k-1]) {
						splice @froms, $k, 1;
						splice @tos,   $k, 1;
						$matches{$name} -= 1;
					}
				}
				#start with true hit length at total coverage then subtract all gaps
				$total_cover  = $pos_ends{$name} - $pos_starts{$name} + 1;

				#subtract  any gap at begining and end using min start and 
				#max end variables from total overage calcuation
				$total_cover -= $s - $pos_starts{$name}; 
				$total_cover -= $pos_ends{$name} - $e;
			
				#if there are still multiple hits after removing total overlaps
				#look for gaps between hits and subtract from total coverage
				for($k = 1; $k <  @froms; $k++) {
					if($tos[$k-1] < $froms[$k]) {
						$total_cover -= $froms[$k] - $tos[$k-1];
					}	
				}			
			}

			printf("%20s	%10d	%10d	%10d	%10d\n", $name, $true_length, $total_cover, $total_over, $matches{$name});
		}
	}
}

