NB. kong project
NB. copyright Lee Render 2009
NB. version 0.9
NB. modules needed
require 'files'
require 'printf'
require 'random'
require 'viewmat'
require 'jfiles'

coclass'psnp'

create=: verb define
NB. member variables ( most get intialised in routines that follow -  don't trust the types )

verbose=: 0
snpdata=:''
snpindices=: ''
seqs=: ''
seqprops=:''
names=:''
data=:''
alignment=:''  NB. unboxed matrix of alignment intermediate file, one row per sequence
props=:''
properties=:''
nice=: ''   NB. final 'nice' summary of data
NUMSNPS=:0
LOGFILE=: 'run.log'
scores=: ''

NB. supplied on command line by client ( snplist.ijs )

INTERMEDIATE_FILE=: ''
GENOTYPES_FILE=: ''
REFERENCE_SEQUENCE_NAME=: ''

)

destroy=: codestroy



NB. getters and setters

get_verbosity =: verb define
	verbose
)



set_verbosity =: verb define
	verbose=: y
)

get_snpindices=: monad define
	snpindices
)
get_scores =: monad define
	scores
)

get_snpdata=: monad define
	snpdata
)

get_data=:  monad define
	NB. snp data for each sequence plus property as last column
	data
)

get_seqprops=: monad define
	seqprops
)

get_names=: monad define
	names
)


NB. useful abbreviations
or=: +.
and=: *.
display =: (1!:2) & 2
le =:  <:

NB. helper functions

NB. logging
log =: 3 : 0
	y fappends <LOGFILE
)

logf =: 4 : 0
	s=. x sprintf y
	log s
	msg s
)


msg=: monad define
if. verbose do.
'%j' printf <y
end.
)

NB. read x delimited string y into a boxed list
parse=: 4 : 0
	( < ;. _2) y,x
)

boxlines =: LF&parse

boxparts =: ' '&parse
genotypeparts =: ';'&parse

count=: dyad define
	val=.   +/  x=y
)

getname=: verb define
	>0{boxparts y
)

getseq =: monad define
	try.
		>1{boxparts y
	catch.
		'ERROR'
	end.
)

getseqname =: monad define
	>0{genotypeparts y
)


getgenotype =: monad define
	>1{genotypeparts y
)

NB. apply counting function to each cell
NB. left argument is the two element vector GFACTOR DFACTOR,
NB. where GFACTOR is the proportion of hyphens we tolerate (e.g. 0.2 ? )
NB. and DFACTOR is the amount by which one letter should dominate the other letters ( e.g. 0.95

get_snps =: dyad define
	'GFACTOR DFACTOR'=. x
	NB. y is the column vector in the aligment (unboxed)
	NB. return 1 if this column has a high proportion of one letter AND doesn't have too high a proportion
	NB. of hyphens
	num =.  #y
	numA=. 'A' count y
	numC=. 'C' count y
	numT=. 'T' count y
	numG=. 'G' count y
	numGarbage =. '-' count y
	garbage =.  numGarbage  % num   NB. proportion of hyphens
	result=. <'NO'

	if. garbage < GFACTOR do.
   		nums=.  numA,numC,numT,numG
   		vec=. nums % num

	 	display vec
	 	bv =.   DFACTOR le vec
	 	display bv
	 	domination=.  1 =   +/ bv    NB. one value dominates

	 	NB. display domination
	 	NB. 'domination value = %j\n' printf domination
	 	if. domination do.
    			NB. value dominates but there should be some other values that are the snp
    			NB. i.e.
    			numgz=. nums > 0
    			spread=. +/ numgz
    			hasmut =.  spread = 2   NB. i.e. only two possible values
    			if. hasmut do.
         			display 'SNP!'
         			NB. display nums
         			NB. get indices of non-zero columns
         			indices =. I. numgz
         			NB. get corresponding letters
         			letters=.  indices { 'ACTG'
         			NB. the numbers of these
         			numsforletters=. indices { nums
         			NB. grade down - Most common codon listed first, SNP second
         			graded=. \: numsforletters 
         			canonical=. graded { letters
         			NB. display canonical
         			result=. <canonical
    			end.

	 	end.

	end.
	result
)

get_canonical_snp =: verb define
	snp_index=. y
	letters=. 'ACTG'
	col=. snp_index { snpindices
	numA=. +/ 'A' = col
	numC=. +/ 'C' = col
	numG=. +/ 'G' = col
	numT=. +/ 'T' = col
	counts=. numA,numC,numG,numT
	graded=. \: counts
	NB. return top two, hence for a SNP: major, minor
	(i.2) { graded { letters	
)


display_snp_values=: verb define
	NB. for each SNP index , show the major minor values
	for_k. snpindices do. 
		values=. get_canonical_snp k
		'SNP;column;%j;%j' logf k;<values
	end.
)	


hassnp =: dyad define
	seq=. > x { seqs
	pos=. y { snpindices
	snp=. > y { snpdata
	mut=. 1 { snp
	codon=. pos { seq
	result =. codon = mut
) 


display_snps_for_seq=: verb define 
		NB. pass in sequence index and log 1-based positions in sequence of any snps
		seqi=. y
		NB. get the corresponding row of the snp table
		snpbv=. seqi { snptable
		snpcols=. I. snpbv
		if. (#snpcols) > 0 do.
			refs=.   seqi col2ref snpcols   NB. the actual positions along the sequence
			name=. > y { names
			snps=.  ": refs
			'sequence;%j;snps;%j' logf name;snps
		end.
)

parsegendata =: monad define
	try.

		seqname=. getseqname y
		prop=. getgenotype y
		props=: props,<prop
		bv=. names = <seqname
		seqindex =. I. bv
		seqprops=: (<prop) (seqindex}) seqprops
	catch.
		NB. if seq not found
		NB. pass
		'error in parsegendata with y = %j' logf <y

	end.

)

pg =: parsegendata @ >
	
snpprop=: dyad define
	NB. count how many of snp x has prop y
	NB. prop column index is the same as NUMSNPS
	NB. ( the column after the last SNP )
	propindex=. NUMSNPS
	snpx=.  (I. >"0 x { |: data ) { data
	propy=. +/ >"0 y&= propindex { |: snpx
)

snpnotprop=: dyad define
	NB. how many of snp x has a prop _other_ than y
	propindex=. NUMSNPS
	snpx=. (I. >"0 x { |: data ) { data
	propnoty=. +/ >"0 y&~: propindex { |: data
)

tanimoto=: dyad define
	NB. dist =: sv1 tanimoto sv2   
	NB. metric to use for two equal sized snp bit-vectors
	or=. -. @: +:   NB. not nor
	intersection=.  +/ x * y
	union=. +/ x or y
	intersection % union
)

D=: monad define
	NB. Simpsons index of diversity
	N=. +/ y
	1 - +/ *: ( y % N)
)

getfullseq=: verb define 
	NB. given sequence index , return string of the sequence
	codons=.   I. 0= '-'= >y{seqs
	codons { >  y { seqs
)

getseqbyname=: verb define
	ind=. I. (<y) = names
	getfullseq ind
)


	

evalsnpprop =: dyad define
	val=. y{properties
	with=. x snpprop val
	without=. x snpnotprop val
	good=: good, y;val;with;without
	vec=. with,without
	graded=. \: vec
	pos=. 0 1
	if. 2=+\graded=pos do.
  		NB. with is bigger or same
  		factor=. 1
	else.
  	factor=. _1
	end.
	g=.factor *  D vec
	'%j evalsnprop %j  = %j     with=%j  without=%j' printf y;x;g;with;without
)	
	

snppropind=: dyad define
	val=. y {properties
	x snpprop val
)

NB. m produces a vector of the property counts of a given snp 
NB. e.g.    m 20 is the list of counts of property values of sequences having snp 20: 0 0 10 .. etc 

m =: 3 : 'y snppropind"0 i.propcount'	

myeval=: dyad  define
NB.   score =: prop_index myeval counts
NB.   counts is a vector of counts of of number of sequences having each snp that ( snp 0, snp 1 .. )  have the property prop_index  
NB.   Essentially take the Euclidean distance to the ideal case
NB.   the closer to  ideal= 0 0 0 0 . . T . 0 0 0   where counti = T i = prop_index
NB.   the better
NB.   sum of counts is 0 we return a bogus high number

t=. +/y
if. t = 0 do.
	'myeval %j  sum is zero!' printf t

	result=. 200
else.
	num=. #y
	NB. "ideal" would be all of them at 
	p=. t ( x } ) ( num # 0 )
	pnorm=. p % t
	ynorm=. y % t
	result=.  %: +/   *: (pnorm - ynorm )
end.
)



widen=: 4 : 0
NB.   score;shortest_index_list =  data widen initial_snp_vector;prop_index
data=. x NB. main boxed data table
sv=. > 0 { y NB. SNP vector of length NUMSNPS
prop_index=. > 1 { y NB. property index of the property we're looking for a minimal test set for
longestlist=. I. sv NB. index list built from SNP vector
TOOBIG=. 5  NB. don't look at more than 5 tests
THRESHOLD=. 0.2
N=. 0
BOGUS=._1
BIGSCORE=. 10000000.00 NB. bogus big score
LEN=.  # longestlist
score=. BIGSCORE	
lowest_score=. score
lowest_il=. _99 _99 NB. bogus initial list

while.  ( N <: LEN ) and  ( N <: TOOBIG ) and  ( score >: THRESHOLD )  do.
	NB. look at one more SNP in the vector
	il=. ( i. N ) { longestlist
	score=. data evalindexlist il;prop_index
	if. score < lowest_score do.
		lowest_score=.score
		lowest_il=. il
	end.
	N=. 1 + N
		
end.

if. lowest_score = BIGSCORE do.
	lowest_score=. BOGUS
end. 

lowest_score;lowest_il
)

indexlistinboxedvec=: 4 : 0
il=.x
vec=.y
*/ >"0 il { vec
)

seqigivenpropi =: 4 : 0
data=.x
val=. y { properties
)

evalindexlist =: 4 : 0
NB. score =: data evalindexlist index_list;property_index
NB. ie work out the diversity index against property value counts  of sequences having the snp vectors
NB. containing the supplied index list il
NB. eg imagine the datatable data was just:
NB.  1 0 1 1 A
NB.  1 1 0 1 A
NB.  0 1 1 0 B
NB.  0 0 1 1 C
NB. one il might be 0 1 3  another would be 1 2. The il 2 3 is "shared" by the first and last rows
data=. x
il=. >0{y
prop_index=. >1{y
prop_value=. prop_index { properties
'rows cols'=. $ data
propcol=. cols - 1
numprops =. #properties

counts=.  numprops # 0

for_r. i.rows do.
	row=. r { data
	if. il indexlistinboxedvec row do.
		boxed_property_value =. propcol { row
		propvec=. boxed_property_value = properties
		prop_ind=.   I. propvec
		oldcount=. prop_ind { counts
		if. prop_ind = prop_index do.
			delta =. 1
		else.
			delta =. 1
		end.

		newcount=. oldcount + delta

		counts =. newcount ( prop_ind } ) counts
		 

	end.
end.
NB. div=. D counts
NB. div
NB. 'prop_index=%j il=%j counts=%j' printf prop_index;il;counts
result=. prop_index myeval counts
NB. 'prop_index=%j il=%j counts=%j score=%j' printf prop_index;il;counts;result
)

evaldata =: dyad define
data=. x
prop_index=. y
'rows cols' =. $ data

best_il=._99 _99
best_score=.10000000  NB. the closer to zero the better ...
for_r. i.rows do.
	sv=. >"0 (i.NUMSNPS) { r { data
	NB. result =. data narrow sv;prop_index
	result =. data widen sv;prop_index
	il=. >1{result
	score=. >0{result
	length=. #il
	if. length >: 1 do.
		if. score < best_score do.
			best_score =. score
			best_il =. il
		end.
	end.
end.
name=. > prop_index { properties
if. ( 0{best_il ) >  0 do.
	snp_cols =. best_il { snpindices
else.
	snp_cols =. _99 _99
end.
NB. 'property;%j;score;%j;snpcols;%j' printf name;best_score;snp_cols
best_score;best_il

)

NB.  position in sequence(1-based) =:  <row of alignment> col2ref col_index ( 0 based ) 
col2ref=: dyad define
	NB.  e.g.  A--CTG col2ref 4 = 3  ( just skip the hyphens )
	NB. XXX what if y{x = '-' ??? 
	y + 1 - +/ '-'=x

)


display_best_lists =: monad define
	scores=. > 0{y
	refseq=. > 1{y
	properties=: >2{y
	refseq=. >3{y
	getpos=. refseq&col2ref
	lists =. 1{ |: scores
	values=.  0{ |: scores 
	positions=. getpos each lists
	numprops=. #properties

	for_i. i.numprops  do.
		propname=. > i { properties
		value=. > i { values
		if. value=0 do.
			list=. >i{ lists
			
			NB. determine the 1-based position of snps in _this_sequence ( not reference , as there isn't any! )
		
			NB. snpcontent=. list { snpdata
			refpositions=. >i{positions
			NB. 'Found a good test for property %j. Test these SNPs ( Positions listed for reference %j): %j' printf propname;refseq;<refpositions
			NB. 'SnpContent: %j' printf <snpcontent	
			cols_string=. ": > list
			ref_pos=. ": refpositions
			'prop;%j;cols;%j;refseq;%j;refpos;%j' logf  propname;cols_string;refseq;ref_pos

		end.
	
	end.
)

run =: monad define
	msg 'Calculating best snp tests for supplied properties - this may take a few minutes ..'
	NB. cut lines into boxes

	INTERMEDIATE_FILE =: > 0{y
	GENOTYPES_FILE =: > 1 { y
	REFERENCE_SEQUENCE_NAME =: > 2 { y
	data =: fread <INTERMEDIATE_FILE
	lines =: boxlines data
	names=: getname each lines
	names =: }: names   NB. drop last (bad record)
	seqs =: getseq each lines
	seqs =: }: seqs
	numseqs =: #seqs
	rows =: # names
	cols =: # >0{ seqs   NB.  should be same length for all
	NB. create matrix
	s =. > ,/  seqs
	alignment=: s

	NB. s should have  shape rows cols - each row is a set if letters ACTG or -
	NB. now we can start to test for SNPs
	NB. box each column by taking transpose
	boxedcols =: <"1 |: s
	NB. look at columns with fewer than 5% hyphens, columns must have 95% or more of one codon
	strict =. 0.05 0.95
	f =: strict&get_snps @ >
	a =: f"0 boxedcols
	notsnp =: <'NO'
	b =: notsnp = a
	snps =: 0 = b
	NB. snpindices = bit vector of columns which have a snp
	snpindices =: I. snps
		
	'snpindices;%j' logf <snpindices


	NB. boxed list of snp data  eg [GT][AC]...   A 
	snpdata =: snpindices { a
	NB. binary table showing which sequences have what snps
	NB. first create iota list for columns and rows
	NB. seqi =. i. ( numseqs)
	NB. snpi =. i. ( # snpindices )
	NB. ensure we're operating on atoms
	NB. r=. 0,0
	NB. f=.  hassnp"r
	NB. produce a table showing what sequences have what snps - must ensure we're operating on atoms, hence then rank modifier
	NB. snptable =. seqi f/ snpi
	snptable =: (i.numseqs) (hassnp"(0 0)) / (i. # snpindices) 
	NB. number of snps = number of colummns
	NUMSNPS =: 1 { $ snptable
	NB. load properties file - ;-delimited
	gendata =: }: boxlines fread <GENOTYPES_FILE
	props =: 0 # <''   NB. empty list
	seqprops=:  numseqs # <'undef'
	numprops =: # props
	pg"0 gendata
	genotypes =: seqprops
	properties=: ~. seqprops

	properties=: (I. (<'undef') ~: properties)  { properties
	propcount =. # properties
	snps =: <"0 snptable
	NB. each row contains snp vector and property value as a boxed list
	data =: snps,.genotypes
	eval =: data&evaldata"0
	msg 'calculating scores - this may take a few minutes ...'
	scores =: eval (i.#properties)
	msg 'finished calculating'


	if. (<REFERENCE_SEQUENCE_NAME) = <'longest' do.	
		NB. use the longest sequence as the reference sequence
		iscodon=. '-'&( -. @ = )
		numcodons=. +/  @ iscodon
		refseqind=.  > 0 { \: numcodons each seqs
	else.
		NB. user supplied reference sequence name
		if. (<REFERENCE_SEQUENCE_NAME) e. names do.
			refseqind=. I. (<REFERENCE_SEQUENCE_NAME) = names
			
		else.
			NB. can't find reference sequence in names list
			msg 'supplied reference sequence name not in aligment'
			exit ''
		end.
	end.


	refseq=.  > refseqind { seqs 
	refseqname=. > refseqind { names
	'reference sequence name  = %j' logf <refseqname

	display_best_lists scores;refseq;properties;refseqname
	display_snps_for_seq"0 i.(#seqs)
	display_snp_values 0
	msg 'done!'
	
	
)


	















 
