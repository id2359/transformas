Instructions

To install:

1. Install the J language interpreter  ( double click on the j602a_win.exe file
supplied.) This installs the required program jconsole.exe.


2. Install Python ( double click on the ActivePython-2.6.4.8-win32-x86.msi file)
This installs the required program python.exe. Make sure that python is
on your path.


3. Edit snplist.cmd and enter the correct path for JROOT
( Which is the folder containing jconsole.exe )


To run:

The program accepts a clustal alignment file and a "properties file"
which simply a text file listing the property / genotype for each
sequence in the alignment ( and no others! ):

E.g

SequenceA;property1
SequenceB;property2
...


Every sequence in the clustal alignment must appear in the list.


To run the program, enter the following at the command prompt:

snplist <clustalalignmentfilename>  <propertyfilename>  <referencesequencename>

OR

snplist <clustalalignmentfilename>  <propertyfilename> 


( which will use the longest sequence in the alignment as a reference.)


The program will create a log file called "run.log"

It _should_ look something like this:

----------------
snpindices;304 306 319 321 323 324 335 343 355 .. < data omitted >
reference sequence name  = BCU70431
prop;5;cols;25;refseq;BCU70431;refpos;26
prop;1;cols;146;refseq;BCU70431;refpos;147
prop;3*1;cols;93;refseq;BCU70431;refpos;94
prop;9;cols;17 36;refseq;BCU70431;refpos;18 37
prop;12;cols;72;refseq;BCU70431;refpos;73
prop;3;cols;63 68;refseq;BCU70431;refpos;64 69
prop;9*1;cols;80 86;refseq;BCU70431;refpos;81 87
prop;thai;cols;27 35;refseq;BCU70431;refpos;28 36
prop;17;cols;40 54;refseq;BCU70431;refpos;41 55
prop;7;cols;41;refseq;BCU70431;refpos;42
prop;11;cols;87;refseq;BCU70431;refpos;88
prop;6;cols;8;refseq;BCU70431;refpos;9
prop;7*g;cols;65 89;refseq;BCU70431;refpos;66 90
prop;16;cols;37;refseq;BCU70431;refpos;38
prop;8;cols;5 95;refseq;BCU70431;refpos;6 96
prop;2;cols;52;refseq;BCU70431;refpos;53
prop;8*1;cols;68 95;refseq;BCU70431;refpos;69 96
prop;13;cols;74;refseq;BCU70431;refpos;75
prop;phym;cols;0;refseq;BCU70431;refpos;1
sequence;AF143795;snps;25
sequence;AY036066;snps;60 121
sequence;cp000614;snps;60 121
sequence;EU563933;snps;26 60 121 127 144
sequence;AF143788;snps;45 52 99 143
sequence;AF143789;snps;35 99
sequence;AF143784;snps;56
sequence;AF143786;snps;45
sequence;AF143787;snps;45 52 54 143
sequence;AF143780;snps;94
sequence;cp000086;snps;28 36 58 59 70 71 78 118
.
.
.


---------------

Explanation:


"snpindices" are the columns  ( starting at column 0 ) of the "intermediate file" ( always called "alignment.im" )
that the program thinks contain a "snp"  ( i.e. mostly one codon but occassionally another one.)

"reference sequence name" is the supplied reference sequence name, or the name of the longest sequence.
"prop" lines show the most informative snps to look for to predict that property.

E.g. in the example data above, the line:


prop;7*g;cols;65 89;refseq;BCU70431;refpos;66 90


means that 

property '7*g' is best tested by looking for snps indicated by the columns 65 and 89 of alignment.im. These are
positions 66 and 90 ( 1-based of sequence BCU70431 ( the reference sequence .)


"sequence" lines list the location ( 1-based ) of each snp in each sequence actually having snps.


E.g. sequence EU563933 has 5 snps .




