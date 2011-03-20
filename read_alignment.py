import os
import sys

line_no = 0
seq_no = 0
datafile  = sys.argv[1]
#datafile = '/home/lee/kong.aln'
aln = {}
names = {}
f = open(datafile,'r')
heading = f.readline().strip()
line_no += 1
heading_parts = heading.split()
num_seqs = int(heading_parts[0])
len_aln =  int(heading_parts[1])
line = f.readline()

while line:
	line = line.strip()
	if line == '':
		seq_no = 0
	else:
		if line_no <= num_seqs + 1:
			# sequence name is first piece
			pieces = line.split()
			name = pieces[0]
			#print name
			alignment = ''.join(pieces[1:])
			#print alignment
			if not aln.has_key(name):
				aln[name] = alignment
				names[seq_no] = name
			else:
				raise Exception("sequence listed twice?")
			
	

		else:
			# we've read the first 'stanza', so pieces are just data with no names
			pieces = line.split()
			alignment = ''.join(pieces)
			# append to existing data
			name = names[seq_no]
			aln[name] = aln[name] + alignment

		seq_no += 1
	

	line = f.readline()
	line_no += 1



f.close()

stream = open("alignment.im","w")

for key in aln.keys():
	l = len(aln[key])
	if l != len_aln:
		print "%s has wrong length %s" % str(l)
	else:
		line = "%s %s\n" % (key,aln[key])
		stream.write(line)


stream.close() 



			
 
			
		


#for line in f.readlines():
#	print line
