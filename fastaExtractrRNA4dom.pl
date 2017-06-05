#!/usr/bin/perl


# this script takes a fasta file with contigs and the output of Barrnap runs for different kingdoms to produce a fasta file with rRNA sequences and a .tab file in the style of prodigal
# 6 inputs: - the fasta file of the contigs
#           - the .gff files from Barrnap (for eukaryota, bacteria, archaea and mitochondria)
#           - the prefix for the output
# 2 outputs: - a fasta file with the rRNA gene sequences, naming is the contig name appended with _r and a continuous number per contig
#	     - a table with contig name rRNA gene name (see above), the sense, length, start position, end position, completeness and kind (16S, 5S etc)

# Anna Heintz-Buschart, June 2016, based on code written in November 2014



use strict;
use warnings;
use Bio::DB::Fasta;
use List::MoreUtils qw(uniq);
use List::Util qw(min max);
use Getopt::Long;

my ($fastaFile,$eukFile,$arcFile,$bacFile,$mitoFile,$output);
GetOptions('f=s' => \$fastaFile,
	   'e=s' => \$eukFile,
	   'a=s' => \$arcFile,
           'b=s' => \$bacFile,
           'm=s' => \$mitoFile,
           'o=s' => \$output);

my $listFile = $output.".tab";
my $geneFile = $output.".fa";


my %rRNAs = (); #0: start, 1: end, 2: sense, 3:kind; 4:eval

open (IN, $bacFile);
while (my $line = <IN>){
    unless ($line =~ /^#/) { 
        chomp $line;
        my($contig, undef, undef, $start, $end, $eval, $sense, undef, $attribute) = split "\t", $line;
        my $completeness = "";
	if ($attribute =~ /partial/){
            $completeness = " incomplete"
        } else {
            $completeness = " complete"
        }
        my @attributes = split ";", $attribute,2;
        my $kind = $attributes[0];
        substr($kind, 0 , 5) = "";
        $kind .= " bacterial".$completeness;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
        push @{$rRNAs{$contig}[2]}, $sense;
        push @{$rRNAs{$contig}[3]}, $kind;
        push @{$rRNAs{$contig}[4]}, $eval;
    }
}
close(IN);

open (IN, $arcFile);
while (my $line = <IN>){
    unless ($line =~ /^#/) { 
        chomp $line;
        my($contig, undef, undef, $start, $end, $eval, $sense, undef, $attribute) = split "\t", $line;
        my $completeness = "";
	if ($attribute =~ /partial/){
            $completeness = " incomplete"
        } else {
            $completeness = " complete"
            }
        my @attributes = split ";", $attribute,2;
        my $kind = $attributes[0];
        substr($kind, 0 , 5) = "";
        $kind .= " archaeal".$completeness;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
        push @{$rRNAs{$contig}[2]}, $sense;
        push @{$rRNAs{$contig}[3]}, $kind;
        push @{$rRNAs{$contig}[4]}, $eval;
    }
}
close(IN);

open (IN, $eukFile);
while (my $line = <IN>){
    unless ($line =~ /^#/) { 
        chomp $line;
        my($contig, undef, undef, $start, $end, $eval, $sense, undef, $attribute) = split "\t", $line;
        my $completeness = "";
	if ($attribute =~ /partial/){
            $completeness = " incomplete"
        } else {
            $completeness = " complete"
            }
        my @attributes = split ";", $attribute,2;
        my $kind = $attributes[0];
        substr($kind, 0 , 5) = "";
        $kind .= " eukaryotic".$completeness;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
        push @{$rRNAs{$contig}[2]}, $sense;
        push @{$rRNAs{$contig}[3]}, $kind;
        push @{$rRNAs{$contig}[4]}, $eval;
    }
}
close(IN);

open (IN, $mitoFile);
while (my $line = <IN>){
    unless ($line =~ /^#/) { 
        chomp $line;
        my($contig, undef, undef, $start, $end, $eval, $sense, undef, $attribute) = split "\t", $line;
        my $completeness = "";
	if ($attribute =~ /partial/){
            $completeness = " incomplete"
        } else {
            $completeness = " complete"
        }
        my @attributes = split ";", $attribute,2;
        my $kind = $attributes[0];
        substr($kind, 0 , 5) = "";
        $kind .= " mitochondrial".$completeness;
        push @{$rRNAs{$contig}[0]}, $start;
        push @{$rRNAs{$contig}[1]}, $end;
        push @{$rRNAs{$contig}[2]}, $sense;
        push @{$rRNAs{$contig}[3]}, $kind;
        push @{$rRNAs{$contig}[4]}, $eval;
    }
}
close(IN);

my %allContigs = ();
my $db = Bio::DB::Fasta->new( $fastaFile );
open(GEN, ">", $geneFile) or die "cannot open $geneFile \n";
open(TAB, ">", $listFile) or die "cannot open $listFile \n";
print TAB "contig\tgene\tsense\tlength\tstart\tend\teval\tkind\n";
my @ids = keys %rRNAs;
foreach my $contig (@ids){
    my @useStarts = ();
    my @useEnds = ();
    my @useSenses = ();
    my @useKinds = ();
    my @useEvals = ();
    if ($#{$rRNAs{$contig}[0]} > 0) {
        my @starts = @{$rRNAs{$contig}[0]};
        my @ends = @{$rRNAs{$contig}[1]};
        my @senses = @{$rRNAs{$contig}[2]};
        my @kinds = @{$rRNAs{$contig}[3]};
        my @evals = @{$rRNAs{$contig}[4]};
        my @idx = sort { $starts[$a]+0 <=> 0+$starts[$b] } 0 .. $#starts;
        @starts = @starts[@idx];
        @ends = @ends[@idx];
        @senses = @senses[@idx];
        @kinds = @kinds[@idx];
        @evals = @evals[@idx];
        my $currentStart = $starts[0];
        my $currentEnd = $ends[0];
        my @currentAllStarts = ();
        my @currentAllEnds = ();
        my @currentSenses = ();
        my @currentKinds = ();
        my @currentEvals = ();
        for (my $i = 0; $i <= $#starts; $i++) { 
            if ($starts[$i] <= $currentEnd){
                $currentEnd = $ends[$i] if ($ends[$i] >= $currentEnd);
                push @currentAllStarts, $starts[$i];
                push @currentAllEnds, $ends[$i];
                push @currentSenses, $senses[$i];
                push @currentKinds, $kinds[$i];
                push @currentEvals, $evals[$i];
            } else {
                my @ind = ();
                my $testEval = $currentEvals[0];
                for (my $m = 0; $m <= $#currentEvals; $m++){
                    if ($testEval == $currentEvals[$m]) {
                        push @ind, $m;
                    } elsif ($testEval > $currentEvals[$m]) {
                        @ind = ($m);
                        $testEval = $currentEvals[$m];
                    }
                }
                unshift @useEvals, $testEval;
                if ($#ind > 0) {
                    my @uniqueSenses = uniq @currentSenses;
                    my @uniqueKinds = uniq @currentKinds;
                    unshift @useSenses, join('',@uniqueSenses);
                    unshift @useKinds, join(';',@uniqueKinds);
                    unshift @useStarts, min @currentAllStarts[@ind];
                    unshift @useEnds, max @currentAllEnds[@ind];
                } else {
                    unshift @useSenses, $currentSenses[$ind[0]];
                    unshift @useKinds, $currentKinds[$ind[0]];
                    unshift @useStarts, $currentAllStarts[$ind[0]];
                    unshift @useEnds, $currentAllEnds[$ind[0]];
                }
                $currentStart = $starts[$i];
                $currentEnd = $ends[$i];
                @currentAllStarts = ($starts[$i]);
                @currentAllEnds = ($ends[$i]);
                @currentSenses = ($senses[$i]);
                @currentKinds = ($kinds[$i]);
                @currentEvals = ($evals[$i]);
            }
        }
        my @ind = ();
        my $testEval = $currentEvals[0];
        for (my $m = 0; $m <= $#currentEvals; $m++){
            if ($testEval == $currentEvals[$m]) {
                push @ind, $m;
            } elsif ($testEval > $currentEvals[$m]) {
                @ind = ($m);
                $testEval = $currentEvals[$m];
            }
        }
        unshift @useEvals, $testEval;
        if ($#ind > 0) {
            my @uniqueSenses = uniq @currentSenses;
            my @uniqueKinds = uniq @currentKinds;
            unshift @useSenses, join('',@uniqueSenses);
            unshift @useKinds, join(';',@uniqueKinds);
            unshift @useStarts, min @currentAllStarts[@ind];
            unshift @useEnds, max @currentAllEnds[@ind];
        } else {
            unshift @useSenses, $currentSenses[$ind[0]];
            unshift @useKinds, $currentKinds[$ind[0]];
            unshift @useStarts, $currentAllStarts[$ind[0]];
            unshift @useEnds, $currentAllEnds[$ind[0]];
        }
    } else {
        push @useStarts, $rRNAs{$contig}[0][0];
        push @useEnds, $rRNAs{$contig}[1][0];
        push @useSenses, $rRNAs{$contig}[2][0];
        push @useKinds, $rRNAs{$contig}[3][0];
        push @useEvals, $rRNAs{$contig}[4][0];
    }
    for (my $j = 0 ; $j <= $#useStarts ; $j++ ) {
        if ($useSenses[$j] =~ /\+-|-\+/) {
            print STDERR "Ambiguous sense on $contig between $useStarts[$j] and $useEnds[$j]. \n";
            next;
        }
        if (exists $allContigs{$contig}) {
            $allContigs{$contig}++
        } else {
            $allContigs{$contig} = 1
        }
        my $gene = join("", $contig, "_r", $allContigs{$contig});
        my $length = 1 + $useEnds[$j] - $useStarts[$j];
        my $sequence = "";
        if ($useSenses[$j] eq "+") {
            $sequence = $db->seq($contig, $useStarts[$j], $useEnds[$j]);
        } else {
            $sequence = $db->seq($contig, $useEnds[$j], $useStarts[$j]);
        }
        if  (!defined( $sequence )) {
            print STDERR "Sequence $contig not found. \n";
            next;
        } else {
            print GEN ">$gene\n", "$sequence\n";
        }
        print TAB join("\t",$contig,$gene,$useSenses[$j],$length,$useStarts[$j],$useEnds[$j],$useEvals[$j],$useKinds[$j]), "\n";

    }
}
close(GEN);
close(TAB);

