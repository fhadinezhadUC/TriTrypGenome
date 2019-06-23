#!/usr/bin/perl -w
# Copyright (C) 2019 David H. Ardell
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# modified from ID_KLD_eps2table on June 18, 2019

use Set::Scalar;
use Getopt::Std;
use vars qw($VERSION $DESC $PROGRAM_NAME); 
$VERSION = do { my @r = (q$Revision: 1.4 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
$DESC    = "Generate text-table input to all.bubble.R from collections of function, ID and KLD logos.\n";
$PROGRAM_NAME   = $0;
$PROGRAM_NAME   =~ s|^\/.*\/||;

$opt_s = 'ACGT'; # The sets of states to include in the table output. To include the gap logo, use '-s "ACGT-"'
$opt_l = undef;  # Increment site-numbers from eps input files by 1 (for LOGOFUN/LOGOTAX input) before mapping to Sprinzl coordinates
&getopts('s:l');

# Format of coordinate-mapping-file
# comma-separated-values
# four coloumns: <site-number-in-eps-file-input>,<x-coord-for-plotting>,<y-coord-for-plotting>,<sprinzl-coordinate>
# an

die "Usage: $PROGRAM_NAME <coordinate-mapping-file> <parasite-functionlogo-file-prefix> <parasite-foreground-IDlogo-file-prefix> <host-foreground-IDlogo-file-prefix> <parasite-foreground-KLDlogo-file-prefix>" unless (@ARGV == 5);

# logo-files should be named with the file-prefixes followed by underscore, state-character, ".eps" as in: <parasite-functionlogo-file-prefix>_A.eps
# Sites in the logo files are assumed to be integers
# To process logos produced from LOGOFUN AND LOGOTAX (which assume that site-labels should start at 0) use option -l. This will add 1 to every site number in the eps-file input.

($coord_mapping_file,$flogo_prefix,$idglogo_prefix,$idllogo_prefix,$kldlogo_prefix) = @ARGV;
@states = split //,$opt_s;
@names = qw/aa coord state fbits fht gainbits gainfht lossbits lossfht convbits convfht x y sprinzl/;
  
$maxcoord = -1000; # THE MAXIMUM SITE-LABEL

open (SKEL,$coord_mapping_file) or die "Could not open $coord_mapping_file\n";
while (<SKEL>){
  chomp;
  my @F = split /,/,$_;
  $X{$F[0]} = $F[1];
  $Y{$F[0]} = $F[2];
  $S{$F[0]} = $F[3];
}
close SKEL;

foreach my $prefix ($flogo_prefix,$idglogo_prefix,$idllogo_prefix,$kldlogo_prefix) {
  foreach my $state (@states) {
    $file = join '', $prefix,"_",$state,".eps";
    open (FILE,$file) or die "Cannot open file $file\n";
    while (<FILE>) {
      if (/^\/barbits -?([\d\.]+)/) {
	$maxbits = $1;   # THE MAXIMUM BITS OF THE LOGO GRAPH
      }
      elsif (/^\/barheight +([\d\.]+)/) {
	$maxheight = $1; # THE HEIGHT IN CENTIMETERS OF THE MAXIMUM BITS
      }
      elsif (/numbering \{\((\d+)\)/) {
	$coord = $1;     # THE GRAPH X-COORDINATE
	$opt_l and $coord++;
	if ($coord > $maxcoord) {
	  $maxcoord = $coord;
	}
      }
      elsif (/^ ([\d\.]+) \((\w)\)/) {
	($ht,$aa) = ($1,$2); # HEIGHT OF LETTER AND ITS SYMBOL
	$totalht += $ht;     # TOTAL HEIGHT OF STACK CALCULATED FROM SYMBOL HEIGHTS
	push @{ $aas{$prefix}{$coord}{$state} },$aa; ## STORES BOTH WHICH CLASSES AND IMPLICITLY, THEIR ORDER
	$ht{$prefix}{$aa} = $ht; # HEIGHT OF STATE
      }
      elsif (/\%\%EndProlog/) {
	$datastart = 1;
      }
      elsif ($datastart and /^grestore$/){ #END OF STACK
	$bits = ($totalht / $maxheight) * $maxbits; # TOTAL BITS IN STACK  
	foreach my $aa (@{ $aas{$prefix}{$coord}{$state} }) { # FOR EACH AA IN STACK
	  $bits{$prefix}{$aa}{$coord}{$state} = $bits; 
	  $fht{$prefix}{$aa}{$coord}{$state} = $ht{$prefix}{$aa} / $totalht; # FRACTIONAL HEIGHT OF LETTER
	  $totalbits{$prefix}{$aa} += $bits;
	}
	$aaset{$prefix}{$coord}{$state} = Set::Scalar->new(@{ $aas{$prefix}{$coord}{$state} }); # SET OF LETTERS IN STACK (COULD USE THRESHHOLD) 
	$totalht = 0;
	$ht = ();
      }
    }
    close FILE;
  }
}


print join "\t",@names,"\n"; ## PRINT HEADER
foreach my $aa (sort keys %{ $totalbits{$flogo_prefix} }) {
  foreach my $coord (1..$maxcoord) {
    foreach my $state (@states) { 
      foreach my $prefix ($flogo_prefix,$idglogo_prefix,$idllogo_prefix,$kldlogo_prefix) {
	unless ($aaset{$prefix}{$coord}{$state} and $aaset{$prefix}{$coord}{$state}->contains($aa)) {
	  $bits{$prefix}{$aa}{$coord}{$state} = $fht{$prefix}{$aa}{$coord}{$state} = 0;
	}
      }
      my $kldbits = ($bits{$kldlogo_prefix}{$aa}{$coord}{$state}); 
      my $kldfht  = ($fht{$kldlogo_prefix}{$aa}{$coord}{$state}); 
      print join "\t",$aa,$coord,$state,$bits{$flogo_prefix}{$aa}{$coord}{$state},$fht{$flogo_prefix}{$aa}{$coord}{$state},$bits{$idglogo_prefix}{$aa}{$coord}{$state},$fht{$idglogo_prefix}{$aa}{$coord}{$state},$bits{$idllogo_prefix}{$aa}{$coord}{$state},$fht{$idllogo_prefix}{$aa}{$coord}{$state},$kldbits,$kldfht,$X{$coord},$Y{$coord},$S{$coord},"\n"
    }
  }
}

