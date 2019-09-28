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


use Getopt::Std;
use vars qw($VERSION $DESC $PROGRAM_NAME); 
$VERSION = do { my @r = (q$Revision: 1.4 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
$DESC    = "\n";
$PROGRAM_NAME   = $0;
$PROGRAM_NAME   =~ s|^\/.*\/||;


 
$LOGO_HEADER = <<'END_HEADER';
%%Pages: atend
%%DocumentFonts:
%%EndComments
/llx  56.7 def
/lly 510.2 def
/urx 896.9 def
/ury 1190.6 def
% * position, samples, information, variance 
%
% logo from 0 to 73

/cmfactor 72 2.54 div def % defines points -> cm conversion
/cm {cmfactor mul} bind def % defines centimeters

% user defined parameters
/lowest 0 def
/highest 73 def
/bar 0 def
/xcorner  4.00000 cm def
/ycorner 26.00000 cm def
/rotation  0.00000 def % degrees
/charwidth  0.36000 cm def
/charwidth2m charwidth 2 mul def
/barheight  5.00000 cm def
/barwidth  0.10000 cm def
/barbits -4.18000 def % bits
/Ibeamfraction  1.00000 def
/barends (b) def
/subticsBig 2 def % sub-tic interval size (1/bits)
/subticsSmall 10 def % sub-tic interval size (1/bits)
/showingbox (n) def
/outline false def
/caps true def
/stacksperline 74 def
/linesperpage 1 def
/linemove  3.00000 def
/numbering true def
/shrinking false def
/edgecontrol (n) def
/edgeleft  2.00000 def
/edgeright  1.00000 def
/edgelow  8.00000 def
/edgehigh  1.00000 def
/shrink  1.00000 def
/ShowEnds (-) def % d: DNA, p: PROTEIN, -: none
/HalfWhiteIbeam false def

/knirhs 1 shrink sub 2 div def
/charwidth4 charwidth 4 div def
/charwidth2 charwidth 2 div def

/outlinewidth {charwidth 32 div} def
/setthelinewidth {% set the linewidth
  outline
    {outlinewidth setlinewidth}
    {1 setlinewidth}
  ifelse
} def
/toggleoutline { % switch the state of outlineing
pop pop pop pop
/outline outline not def
setthelinewidth
} def

% define fonts
/ffss {findfont fontsize scalefont setfont} def
/FontForStringRegular {/Times-Bold       ffss} def
/FontForStringItalic  {/Times-BoldItalic ffss} def
/FontForLogo          {/Helvetica-Bold   ffss} def
/FontForPrime         {/Symbol           ffss} def
/FontForSymbol        {/Symbol           ffss} def

% Set up the font size for the graphics
/fontsize charwidth def

% movements to place 5' and 3' symbols
/fivemovex {0} def
/fivemovey {(0) charparams lx ux sub 3 mul} def
/threemovex {(0) stringwidth pop 0.5 mul} def
/threemovey {fivemovey} def
/prime {FontForPrime (\242) show FontForStringRegular} def

% make italics possible in titles
/IT {% TRstring ITstring IT -
  exch show
  FontForStringItalic
  show
  FontForStringRegular
} def


% make symbols possible in titles
/SY {% TRstring SYstring SY -
  exch show
  FontForSymbol
  show
  FontForStringRegular
} def

%(*[[ This special comment allows deletion of the repeated
% procedures when several logos are concatenated together
% See the censor program.

/charparams { % char charparams => uy ux ly lx
% takes a single character and returns the coordinates that
% defines the outer bounds of where the ink goes
  gsave
    newpath
    0 0 moveto
    % take the character off the stack and use it here:
    true charpath 
    flattenpath 
    pathbbox % compute bounding box of 1 pt. char => lx ly ux uy
    % the path is here, but toss it away ...
  grestore
  /uy exch def
  /ux exch def
  /ly exch def
  /lx exch def
} bind def

/dashbox { % xsize ysize dashbox -
% draw a dashed box of xsize by ysize (in points)
  /ysize exch def % the y size of the box
  /xsize exch def % the x size of the box
  1 setlinewidth
  gsave
    % Define the width of the dashed lines for boxes:
    newpath
    0 0 moveto
    xsize 0 lineto
    xsize ysize lineto
    0 ysize lineto
    0 0 lineto
    [3] 0 setdash
    stroke
  grestore
  setthelinewidth
} bind def

/boxshow { % xsize ysize char boxshow
% show the character with a box around it, sizes in points
gsave
  /tc exch def % define the character
  /ysize exch def % the y size of the character
  /xsize exch def % the x size of the character
  /xmulfactor 1 def /ymulfactor 1 def

  % if ysize is negative, make everything upside down!
  ysize 0 lt {
    % put ysize normal in this orientation
    /ysize ysize abs def
    xsize ysize translate
    180 rotate
  } if

  shrinking {
    xsize knirhs mul ysize knirhs mul translate
    shrink shrink scale
  } if

  2 {
    gsave
    xmulfactor ymulfactor scale
    tc charparams
    grestore

    ysize % desired size of character in points
    uy ly sub % height of character in points
    dup 0.0 ne {
      div % factor by which to scale up the character
      /ymulfactor exch def
    } % end if
    {pop pop}
    ifelse

    xsize % desired size of character in points
    ux lx sub % width of character in points
    dup 0.0 ne {
      div % factor by which to scale up the character
      /xmulfactor exch def
    } % end if
    {pop pop}
    ifelse
  } repeat

  % Adjust horizontal position if the symbol is an I
  tc (I) eq {charwidth 2 div % half of requested character width
             ux lx sub 1 div % half of the actual character
                sub      0 translate} if
  % Avoid x scaling for I
  tc (I) eq {/xmulfactor 2 def} if

  /xmove xmulfactor lx mul neg def
  /ymove ymulfactor ly mul neg def

  newpath
  xmove ymove moveto
  xmulfactor ymulfactor scale

  outline {  % outline characters:
setthelinewidth
    tc true charpath
    gsave 1 setgray fill grestore
    clip stroke
}
  { % regular characters
    tc show
  }
  ifelse
grestore
} def

/numchar{ % charheight character numchar
% Make a character of given height in cm,
% then move vertically by that amount
  gsave
    /char exch def
    /charheight exch cm def
    /visible true def % most characters are visible
    char (K) eq {0 0 1 setrgbcolor} if
    char (R) eq {0 0 1 setrgbcolor} if
    char (H) eq {0 0 1 setrgbcolor} if
    char (k) eq {0 0 1 setrgbcolor} if
    char (r) eq {0 0 1 setrgbcolor} if
    char (h) eq {0 0 1 setrgbcolor} if
    char (D) eq {1 0 0 setrgbcolor} if
    char (E) eq {1 0 0 setrgbcolor} if
    char (d) eq {1 0 0 setrgbcolor} if
    char (e) eq {1 0 0 setrgbcolor} if
    char (N) eq {  0.6000 0   0.8500 setrgbcolor} if
    char (Q) eq {  0.6000 0   0.8500 setrgbcolor} if
    char (n) eq {  0.6000 0   0.8500 setrgbcolor} if
    char (q) eq {  0.6000 0   0.8500 setrgbcolor} if
    char (F) eq {1 0 1 setrgbcolor} if
    char (Y) eq {1 0 1 setrgbcolor} if
    char (W) eq {1 0 1 setrgbcolor} if
    char (f) eq {1 0 1 setrgbcolor} if
    char (y) eq {1 0 1 setrgbcolor} if
    char (w) eq {1 0 1 setrgbcolor} if
    char (G) eq {0   0.7000 0 setrgbcolor} if
    char (A) eq {0   0.7000 0 setrgbcolor} if
    char (S) eq {0   0.7000 0 setrgbcolor} if
    char (T) eq {0   0.7000 0 setrgbcolor} if
    char (g) eq {0   0.7000 0 setrgbcolor} if
    char (a) eq {0   0.7000 0 setrgbcolor} if
    char (s) eq {0   0.7000 0 setrgbcolor} if
    char (t) eq {0   0.7000 0 setrgbcolor} if
    char (C) eq {1   0.8500 0 setrgbcolor} if
    char (c) eq {1   0.8500 0 setrgbcolor} if
    char (P) eq {0 1 1 setrgbcolor} if
    char (p) eq {0 1 1 setrgbcolor} if
    char (X) eq {0 0 0 setrgbcolor} if
    char (M) eq {0 0 0 setrgbcolor} if
    char (I) eq {0 0 0 setrgbcolor} if
    char (L) eq {0 0 0 setrgbcolor} if
    char (V) eq {0 0 0 setrgbcolor} if
    char (x) eq {0 0 0 setrgbcolor} if
    char (m) eq {0 0 0 setrgbcolor} if
    char (i) eq {0 0 0 setrgbcolor} if
    char (l) eq {0 0 0 setrgbcolor} if
    char (v) eq {0 0 0 setrgbcolor} if
    char (*) eq {1 0 0 setrgbcolor} if
    char (\\) eq {0 1 0 setrgbcolor} if
    char (\%) eq {0 0 1 setrgbcolor} if
    char (+) eq {1 0 1 setrgbcolor} if
     visible {
       % implement boxes, fill and characters:
       showingbox (s) eq
       showingbox (f) eq
       or
       {gsave
           ly lx
           ly charwidth add
           lx charheight add
           boxsymbol
           clip
           showingbox (f) eq
           {fill}  
           {gsave 0 setgray stroke grestore  
            charwidth charheight char boxshow
           }
           ifelse
       grestore
       }
       {charwidth charheight char boxshow}
       ifelse
     } if % visibility control
  grestore
  charheight abs 1 gt {0 charheight abs translate} if
} bind def

/Ibar{
% make a horizontal bar
gsave
  newpath
    charwidth4 neg 0 moveto
    charwidth4 0 lineto
  stroke
grestore
} bind def

/Ibeam{ % height Ibeam
% Make an Ibeam of twice the given height, in cm
  /height exch cm def
  /heightDRAW height Ibeamfraction mul def
  1 setlinewidth
     HalfWhiteIbeam outline not and
     {0.75 setgray} % grey on bottom
     {0 setgray} % black on bottom
  ifelse
  gsave
    charwidth2 height neg translate
    Ibar
    newpath
      0 0 moveto
      0 heightDRAW rlineto
    stroke
    0 setgray % black on top
    newpath
      0 height moveto
      0 height rmoveto
      currentpoint translate
    Ibar
    newpath
      0 0 moveto
      0 heightDRAW neg rlineto
      currentpoint translate
    stroke
  grestore
  setthelinewidth
} bind def

/makenumber { % number makenumber
% make the number
gsave
  shift % shift to the other side of the stack
  90 rotate % rotate so the number fits
  dup stringwidth pop % find the length of the number
  neg % prepare for move
  charwidth (0) charparams uy ly sub % height of numbers
  sub 2 div %
  moveto % move back to provide space
  show
grestore
} bind def

/shift{ % move to the next horizontal position
charwidth 0 translate
} bind def

/bar2 barwidth 2 div def
/bar2n bar2 neg def
/makebar { % make a vertical bar at the current location
gsave
   bar2n 0 moveto
   barwidth 0 rlineto
   0 barheight rlineto
   barwidth neg 0 rlineto
   closepath
   fill
grestore
} def

% definitions for maketic
/str 10 string def % string to hold number
% points of movement between tic marks:
% (abs protects against barbits being negative)
/ticmovement barheight barbits abs div def

/maketic { % make tic marks and numbers
% define tic mark to be the width of the number 4:
(4) stringwidth pop
/ticwidth exch def % width of tic (as a dash) to show
gsave
  % initial increment limit proc for
  0 1 barbits abs cvi
  {/loopnumber exch def

    % convert the number coming from the loop to a string
    % and find its width
    loopnumber 10 str cvrs
    /stringnumber exch def % string representing the number

    stringnumber stringwidth pop
    /numberwidth exch def % width of number to show

    /halfnumberheight
      stringnumber charparams % capture sizes
      uy ly sub 2 div
    def


    numberwidth % move back width of number
    neg loopnumber ticmovement mul % shift on y axis
    halfnumberheight sub % down half the digit

    moveto % move back the width of the string

    ticwidth neg 0 rmoveto % move back the width of the tic

    stringnumber show

    % now show the tic mark
    0 halfnumberheight rmoveto % shift up again
    ticwidth 0 rlineto
    stroke
  } for
grestore

% do additional BIG tic marks.  subtics is user defined
  % initial increment limit proc for
gsave
  0 1 barbits subticsBig mul abs cvi
  {/bitnumber exch subticsBig div subticsBig div def
    0
    neg bitnumber ticmovement mul subticsBig mul % shift on y axis
    moveto
    ticwidth neg 0 rlineto
    stroke
  } for
/subticsBig 2 def % sub-tic interval size (1/bits)
% do additional SMALL tic marks.  subticsSmall is user defined
/ticwidth ticwidth 2 div def % halve the ticwidth
  % initial increment limit proc for
gsave
  0 1 barbits subticsSmall mul abs cvi
  {/bitnumber exch subticsSmall div subticsSmall div def
    0
    neg bitnumber ticmovement mul subticsSmall mul % shift on y axis
    moveto
    ticwidth neg 0 rlineto
    stroke
  } for
grestore
gsave
  /labelstring (bits) def
  numberwidth neg 2.5 mul
  barheight
  labelstring stringwidth pop
  sub 2 div
  translate
  90 rotate
  0 0 moveto
  labelstring show
grestore
} def

/degpercycle 360 def
 
/drawcosine {% amplitude  phase  wavelength  base
%              xmin ymin xmax ymax step
%              dashon dashoff dashoffset thickness
%              cosine -
% draws a cosine wave with the given parameters:
% amplitude (points): height of the wave
% phase (points): starting point of the wave
% wavelength (points): length from crest to crest
% base (points): lowest point of the curve
% xmin ymin xmax ymax (points): region in which to draw
% step steps for drawing a cosine wave
% dashon if greater than zero, size of dashes of the wave (points)
% dashon dashing on interval (points)
% dashoff dashing off interval (points)
% dashoffset offset for dashing (points)
% thickness if greater than zero, thickness of wave (points)
% use dashon and dashoff as blank and dashoffset as 0 for solid line
% See PostScrirt Language Reference Manual 2nd ed p. 500 on dash.

  /thickness exch def
  /dashoffset exch def
  /dashoff exch def
  /dashon exch def
  /step exch def
  /ymax exch def
  /xmax exch def
  /ymin exch def
  /xmin exch def
  /base exch def
  /wavelength exch def
  /phase exch def
  /amplitude exch def
  % fun := amplitude*cos( ((-y-phase)/wavelength)*360) + base
  /fun {phase sub wavelength div degpercycle mul cos
           amplitude mul base add} def

  gsave
    /originallinewidth currentlinewidth def
    thickness 0 gt {thickness setlinewidth} if

    % Force the curve to fit into the region specified:
    newpath
    xmin ymin moveto
    xmax ymin lineto
    xmax ymax lineto
    xmin ymax lineto
    closepath
    clip

    newpath
    xmin dup fun moveto
    % go to xmin-1 and xmax+1 to make sure we overlap the
    % next wave if there is one.  The clip above ensures that it
    % goes no further than requested. 
    % loop from xmin-1 to xmax+1 by step:
    xmin 1 sub step xmax 1 add {dup fun lineto} for
    % turn dash on if dashon is positive
    dashon 0 gt {[dashon cvi dashoff cvi] dashoffset setdash} if
    stroke

    originallinewidth setlinewidth
  grestore
} bind def

/circlesymbol { % x y radius circlesymbol - (path)
newpath 0 360 arc closepath} bind def

/sqrt3 3 sqrt def
/trianglesymbol { % x y radius trianglesymbol - (path)
/r exch def
/sqrt3r sqrt3 r mul def
translate
120 rotate
0 r translate
-120 rotate
newpath
0 0 moveto
sqrt3r 0 lineto
-300 rotate
sqrt3r 0 lineto
closepath} bind def

/squaresymbol { % x y side squaresymbol - (path)
/side exch def
translate
side 2 div neg dup translate
newpath
0 0 moveto
0 side lineto
side side lineto
side 0 lineto
closepath} bind def

/linesymbol { % x1 y1 x2 y2 linesymbol - (path)
/y2 exch def
/x2 exch def
/y1 exch def
/x1 exch def
newpath
x1 y1 moveto
x2 y2 lineto
} bind def

/boxsymbol { % x1 y1 x2 y2 boxsymbol - (path)
/y2 exch def
/x2 exch def
/y1 exch def
/x1 exch def
newpath
x1 y1 moveto
x2 y1 lineto
x2 y2 lineto
x1 y2 lineto
closepath
} bind def

% The following special comment allows deletion of the repeated
% procedures when several logos are concatenated together
% See the censor program.
%]]%*)

/startpage { % start a page
  save % [ startpage
  % set the font used in the title strings
  FontForStringRegular
  gsave % [ startpage
  xcorner ycorner translate
  rotation rotate
  % create the user defined strings
  gsave
    /stringscale  2.00000 def
     0.00000 cm -1.00000 cm moveto
    stringscale stringscale scale
    ()
    show
  grestore
  gsave
    % string number 1
    % center the string
    /stringscale  2.00000 def
    ()
    stringwidth pop
    stringscale mul neg
    stacksperline charwidth mul
    add 2 div
    -1.00000 cm moveto
    stringscale stringscale scale
    ()
    show
  grestore
  % now move up to the top of the top line:
  0 linesperpage linemove barheight mul mul translate

  % set the font used in the logos
  FontForLogo
} def

%(*[[ This special comment allows deletion of the repeated
% procedures when several logos are concatenated together
% See the censor program.

/endpage { % end a page
  grestore % ] endpage
  showpage % REMOVE FOR PACKAGING INTO ANOTHER FIGURE
  restore % ] endpage
} def

/showleftend {
gsave
 charwidth neg 0 translate
 fivemovex fivemovey moveto ShowEnds (d) eq {(5) show prime} if
 ShowEnds (p) eq {(N) show} if
grestore
} def

/showrightend {
gsave
 threemovex threemovey moveto ShowEnds (d) eq {(3) show prime} if
 ShowEnds (p) eq {(C) show} if
grestore
} def

/startline{ % start a line
% move down to the bottom of the line:
  0 linemove barheight mul neg translate
  gsave % [ startline
  % put a bar on the left side:
  barends (b) eq barends (l) eq or {
    maketic % maketic.startline
    gsave
      bar2n 0 translate % makebar.startline
      makebar % makebar.startline
    grestore
  } if
  showleftend
} def

/endline{ % end a line
  showrightend
  % put a bar on the right side:
  barends (b) eq barends (r) eq or {
    gsave
      bar2 0 translate % makebar.endline
      makebar % makebar.endline
    grestore
  } if
  grestore % ] startline
} def

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@ End of procedures @@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% The following special comment allows deletion of the repeated
% procedures when several logos are concatenated together
% See the censor program.
%]]%*)

END_HEADER


die "Usage: $PROGRAM_NAME *.eps" unless (@ARGV > 0);

foreach my $file (@ARGV) {
  $newfile = join '', $file,".orig";
  rename $file,$newfile;
  open (FILE,"<$newfile") or die "Cannot open file $newfile\n";
  open (OUT,">$file") or die "Cannot open file $file\n";

  if ($file =~ /[ACGU-][ACGU-]_/){
    $llx = 100;
    $lly = 695;
    $urx = 338;
    $ury = 885;
  }
  else{
    $llx = 100;
    $lly = 710;
    $urx = 875;
    $ury = 885;
  }
  print OUT <<"TOP";
%!PS-Adobe-2.0 EPSF-2.0
%%Title: makelogo 9.34
%%Creator: Tom Schneider, toms\@ncifcrf.gov
%%BoundingBox:    $llx   $lly   $urx  $ury
TOP
  
  print OUT $LOGO_HEADER;
  print $skip = 1;
  while (<FILE>) {
    if (/%%EndProlog/) {
      undef $skip;
    }
    print OUT $_ unless ($skip);
  }
  close FILE;
  close OUT;
};
1;

