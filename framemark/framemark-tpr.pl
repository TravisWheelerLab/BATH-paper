use strict;
use warnings;
use List::Util qw[min max];
my @hit;
my $i;
my $j;
my $k;

my $usage = "Usage: perl framemark-coverage.pl <merfile> \n";

my $merfile = shift;
my $tool = shift;

my %num_targets = ();
my %tp_found = ();
my %tp_before_fp = ();
my @fields;
my $fam;
my $line;
my $type;

#get the total numner of target sequences (true positives) for each familiy
open(MER, $merfile) || die "ERROR, couldn't open $merfile";
while ($line = <MER>)
{
    chomp $line;
    if    ($line =~ m/^\#/) { next; }
    elsif ($line =~ m/^\=/) { next; }
    else  {
        @fields   = split(' ', $line, 6);
        $fam = $fields[0];
        if ($fam =~ m/^\*/) { next; }
        $num_targets{$fam} = $fields[2]; 
    }
}
close(MER);

# count all the true positives found for each familiy
# Once the first false positive is found record the tp found up to that point
open(MER, $merfile) || die "ERROR, couldn't open $merfile";
while ($line = <MER>)
{
    chomp $line;
    if    ($line =~ m/^\=/) {  
        @fields   = split(' ', $line, 6);
	$fam = $fields[2];
	$type = $fields[5];
	if ($type =~ m/\+/) {
            if (exists $tp_found{$fam}) { $tp_found{$fam}++;   } 
            else                        { $tp_found{$fam} = 1; }
         
            if    (exists $tp_before_fp{$fam})            { next; }   
            elsif ($tp_found{$fam} == $num_targets{$fam}) { $tp_before_fp{$fam} = $tp_found{$fam}; }
        }
        else {
            if    (exists $tp_before_fp{$fam}) { next; }
	    elsif (exists $tp_found{$fam})     { $tp_before_fp{$fam} = $tp_found{$fam}; }
            else                               { $tp_before_fp{$fam} = 0; }
        }
    }
}

close(MER);

printf("%10s %20s %10s %10s %16s\n", "# tool", "family", "total tp", "tp found", "tp before 1st fp");

foreach $fam (keys %num_targets) {
    printf("%10s %20s %10d ", $tool, $fam, $num_targets{$fam});

    if (exists $tp_found{$fam}) { printf("%10d ", $tp_found{$fam}); }
    else                        { printf("%10d ", 0); }

    if (exists $tp_before_fp{$fam}) { printf("%10d\n", $tp_before_fp{$fam}); }
    else                            { printf("%10d\n", 0); }
}

