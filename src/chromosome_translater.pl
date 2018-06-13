#!/usr/bin/perl

# Title: table_cutter.pl
# Author: Jochen Bick
# Description:
# Translates chromosome IDs to numbers

use strict;
#use Data::Dumper;
use Getopt::Long;
my $gff;
my $chr = 0;
my $out = "chr_trans.gff";

GetOptions ("chr=s" => \$chr, # string
	    "gff=s" => \$gff, # string
	    "out=s" => \$out) # string
    or die("Error in command line arguments\n");

my %chr_hash;
#my $comment;
#open(GFF, "<$gff") or die "Cannot open: $!";
#open(CHR, "<$chr") or die "Cannot open: $!";
#open(OUT, ">$out") or die "Cannot open: $!";

my @chrfiles = split(/,/, $chr);
#print Dumper @chrfiles;
if($chr){
    foreach my $chrfile (@chrfiles){
	open(CHR, "<$chrfile") or die "Cannot open: $!";
	while (my $line = <CHR>){
	    chomp($line);
	    next if($line =~ /^#/);
	    my @elements = split(/\t/, $line);
	    if($elements[0] =~ /Un/){
		$elements[3] =~ s/\..*//;
		$chr_hash{$elements[1]} = "chrUn_".$elements[3];
		next;
	    }
	    #$elements[0] =~ s/(\w)\w/$1/; # to change MT to M only
	    $elements[0] =~ s/MT/M/; # to change MT to M only	 
	    $chr_hash{$elements[1]} = "chr".$elements[0];
	}
	close(CHR);
    }
    open(GFF, "<$gff") or die "Cannot open: $!";
    open(OUT, ">$out") or die "Cannot open: $!";
    
    while (my $line = <GFF>){
	chomp($line);
	if($line =~ /^#/){ # for lines with comment
	    print OUT $line."\n";
	    next;
	}
	my @lines = split(/\t/, $line);
	if($chr_hash{$lines[0]}){
	    $lines[0] = $chr_hash{$lines[0]};
	}
	my $geneid = $lines[8];
	my $genbank = $lines[8];
	if($lines[8] !~ /GeneID:(\d+)/){
	    next;
	}
	$geneid =~ s/.*GeneID:(\d+).*/$1/;
	if($genbank =~ /.*Genbank:(\w\w_\d+\.\d+).*/){
	    $genbank = $1;
	}else{
	    $genbank = "-";
	}
	$lines[8] .= ";geneid=$geneid;genbank=$genbank";

	print OUT join("\t", @lines);
	print OUT "\n";
    }
    close(OUT);
    close(GFF);
}else{
    open(GFF, "<$gff") or die "Cannot open: $!";
    open(OUT, ">$out") or die "Cannot open: $!";
    
    while (my $line = <GFF>){
	chomp($line);
	if($line =~ /^#/){ # for lines with comment
	    print OUT $line."\n";
	    next;
	}
	my @lines = split(/\t/, $line);
	if($lines[0] =~ /^[\w\d]{1,2}$/){
	   # $lines[0] =~ s/(\w)\w/$1/; # to change MT to M only	    
	    $lines[0] =~ s/^MT$/M/; # to change MT to M only
	    $lines[0] = "chr".$lines[0];
	}else{
	    $lines[0] =~ s/\..*//;
	    $lines[0] = "chrUn_".$lines[0];
	}
	
	print OUT join("\t", @lines);
	print OUT "\n";
    }
    close(OUT);
    close(GFF);
}

