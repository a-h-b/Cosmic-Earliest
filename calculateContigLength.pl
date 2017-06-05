#!/usr/bin/perl

#contributed by Patrick May


use strict;

my $contigs=$ARGV[0];

my %contigs=();

my ($id);

print "sequenceID\tlength\tGCperc\n"; 

open(CON,$contigs) or die $!;
while(my $str=<CON>){
    chomp($str);
    next if length($str)==0;
    if($str=~/>(\S+)/){
	$id=$1;
	die "$id already exists\n" if exists($contigs{$id});
    }else{
	$contigs{$id}{'seq'}.=$str;
    }
} 
close(CON);


foreach $id (keys(%contigs)){
    next if length($id)==0;
    my $seq=$contigs{$id}{'seq'};
    my $length=length($seq); 
    if ($length eq 0){
	print STDERR "ERROR: $id $length $seq\n";
	next;
    }
    my $gc = $seq =~ tr/GC//;
    my $gcp = sprintf("%.4f",$gc/$length);
    $id=~s/\s+/_/g;
    print "$id\t$length\t$gcp\n"; 
}
