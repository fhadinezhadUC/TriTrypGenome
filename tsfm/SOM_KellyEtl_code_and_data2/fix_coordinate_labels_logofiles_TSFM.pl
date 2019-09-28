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

# modified from fix_coordinate_labels_logofiles.pl on June 26, 2019 to handle Class-Informative Base-Pairs.
# modified from ID_KLD_eps2table on June 18, 2019

# RUN THIS AFTER GENERATING BUBBLEPLOTS FROM LOGOS FOR FINAL PRESENTATION OF LOGOS
# THIS WILL CAUSE THE LOGO INPUT TO BREAK THE BUBBLEPLOT SCRIPT

use Set::Scalar;
use Getopt::Std;
use vars qw($VERSION $DESC $PROGRAM_NAME); 
$VERSION = do { my @r = (q$Revision: 1.4 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
$DESC    = "\n";
$PROGRAM_NAME   = $0;
$PROGRAM_NAME   =~ s|^\/.*\/||;


# Format of coordinate-mapping-file
# comma-separated-values
# four coloumns: <site-number-in-eps-file-input>,<x-coord-for-plotting>,<y-coord-for-plotting>,<sprinzl-coordinate>
# 

die "Usage: $PROGRAM_NAME <coordinate-mapping-file> *.eps" unless (@ARGV > 1);

# logo-files should be named with the file-prefixes followed by underscore, state-character, ".eps" as in: <parasite-functionlogo-file-prefix>_A.eps
# Site coords in the logo files are assumed to be integers starting from 0
# To process logos produced from LOGOFUN AND LOGOTAX (which assume that site-labels should start at 0) use option -l. This will add 1 to every site number in the eps-file input.

# NA,7.250,8.125,-1
# 1,6.875,8.875,1
# 2,6.500,8.500,2
# 3,6.125,8.875,3
# 4,5.750,8.500,4
# 5,5.375,8.875,5
# 6,5.000,8.500,6
# 7,4.625,8.875,7
# 8,4.625,7.375,8
# 9,5.000,7.000,9


$coord_mapping_file = shift @ARGV;
  
open (SKEL,$coord_mapping_file) or die "Could not open $coord_mapping_file\n";
while (<SKEL>){
  chomp;
  next if (/^NA/);
  my @F = split /,/,$_;
  $S{$F[0]} = $F[3];
}
close SKEL;

foreach my $file (@ARGV) {
  $newfile = join '', $file,".orig";
  rename $file,$newfile;
  open (FILE,"<$newfile") or die "Cannot open file $newfile\n";
  open (OUT,">$file") or die "Cannot open file $file\n";
  while (<FILE>) {
    if (/numbering \{\((\d+)\)/) {
      $coord = $1;		# THE GRAPH X-COORDINATE
      $mapped = $coord;
      $mapped++;
      if (exists  $S{$mapped}) {
	my $S = $S{$mapped};
	$_ =~ s/\{\(\d+/{($S/;
      } else {
	print STDERR "warning: eps file contains coordinate label $coord we mapped to $mapped, but not found in L-SKEL file.\n";
      }

    }
    elsif (/numbering \{\(\((\d+), (\d+)\)\)/) {
      $coord1 = $1;
      $coord2 = $2;# THE GRAPH X-COORDINATE
      $mapped1 = $coord1 + 1;
      $mapped2 = $coord2 + 1;
      if (exists  $S{$mapped1} and exists $S{$mapped2}) {
	my $S1 = $S{$mapped1};
	my $S2 = $S{$mapped2};
	$_ =~ s/\{\(\(\d+/{(($S1/;
	$_ =~ s/\d+\)\)/$S2))/;
      } else {
	print STDERR "warning: eps file contains coordinate label $coord1 and $coord2 we mapped to $mapped1 and $mapped2, but not found in L-SKEL file.\n";
      }

    }
    print OUT $_;
  }
  close FILE;
  close OUT;
}


