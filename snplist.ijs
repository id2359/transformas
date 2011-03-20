load 'psnp.ijs'
intermediate_file =: > 2 { ARGV_j_
properties_file =:   > 3 { ARGV_j_
reference_sequence =: > 4 { ARGV_j_
snplister =: conew 'psnp'
create__snplister 0
set_verbosity__snplister 1
main=: verb define
try.
	data =:   run__snplister intermediate_file;properties_file;<reference_sequence
	exit 0
catch.
	exit 1
end.
)

main ''

