# Sorts reads with quaternary Hamming codes and de
# set of valid hamming codes
# sys.argv[1] tab delimited barcode file
# attach degenerate sequence to read id "_NNNN"

import argparse, sys, gzip, itertools, glob, pdb
parser = argparse.ArgumentParser(description='Sort barcodes from gzipped fastq file')
parser.add_argument('barcodes', help='tab-delimited file with columns for library and the two barcodes')
parser.add_argument('fastq', help='input fastq file, or glob')
parser.add_argument('prefix', help='prefix of the output file')

args = parser.parse_args()


def fix(barcode):
	"""
	Accepts AC transposed versions of H4_74 barcodes Table S2 Bystrykh 2012
	Assumes all uppercase sequences with ACTGN alphabet.
	If there is more than one N in the sequence, the function returns False
	"""
	def transpose(seq):
		#AC transpose
		seq = seq.replace('C','X')
		seq = seq.replace('A','C')
		return(seq.replace('X','A'))

	nuc2code = dict(zip(["A","C","G","T","N"],range(4)+[0])) # convert characters to numeric, give N values a value of A
	base=['A','C','G','T'] 
	p = [0]*3
	HQ74 = []
	barcode = list(transpose(barcode))
	#convert code to numbers, and check to see if there is too much erroneous data
	bad_count = 0
	for n in barcode:
		HQ74.append(nuc2code[n])
		if n == "N":
			bad_count += 1
		# check to see if barcode has too much missing data
		if bad_count > 1:
			return(False)

	# calculate the checksums
	p[0]=(HQ74[0]+HQ74[2]+HQ74[4]+HQ74[6])%4
	p[1]=(HQ74[1]+HQ74[2]+HQ74[5]+HQ74[6])%4
	p[2]=(HQ74[3]+HQ74[4]+HQ74[5]+HQ74[6])%4 
	p_max=max(p)
	
	if p_max>0:
		#convert checksums to binary
		for i in range (3):
			if p[i]>0:
				p[i]=1
				
		#compute location of error
		err_pos = 0 | p[2]*4 | p[1]*2 | p[0] - 1

		#wrong value min error type %4 should give the right answer
		corr_value=(HQ74[err_pos]-p_max)%4
		barcode[err_pos] = base[corr_value]
	return(transpose("".join(barcode)))


def parse_fastq(handle):
	rec = []
	count = 0
	for line in handle:
		if line:
			rec.append(line)
			count +=1
			if count == 4:
				yield rec
				rec = []
				count = 0

#read barcodes and open a output files
outfiles = {}
codes = {}
degen = 4
#sample id, barcode 1, barcode 2

outfiles['unsorted'] = gzip.open(args.prefix + '_unsorted.fq.gz',"w")
for line in open(args.barcodes):
	if line.rstrip():
		(name,bc1,bc2) = line.split()
		outfiles[name] = gzip.open(args.prefix + name + ".fq.gz","w")
		codes[(bc1,bc2)] = name

unsorted = 0
total = 0

for f in glob.glob(args.fastq):
	f = gzip.open(f)
	for rec in parse_fastq(f):
		total += 1
		header = rec[0]
		barcode = header.rstrip().split()[-1].split(":")[-1]
		bc1 = barcode[0:7]
		Ns = barcode[7:7+degen]
		# attach degenerate tag to read id
		readId = rec[0].split(" ")
		readId[0] = readId[0] + "_" + Ns
		rec[0] = " ".join(readId)
		bc2 = barcode.split("+")[-1][:7] # assume barcodes are separated by a space
		#try to find the correct barcoce, trying fix
		name = codes.get((bc1,bc2))
		if not name:
			fix_bc1 = fix(bc1)
			#if permanently broken, continue
			if not fix_bc1:
				unsorted += 1
				outfiles['unsorted'].write("".join(rec))
				continue
			name = codes.get((fix_bc1,bc2))
			if not name:
				fix_bc2 = fix(bc2)
				if not fix_bc2:
					unsorted += 1
					outfiles['unsorted'].write("".join(rec))
					continue
				name = codes.get((bc1,fix_bc2))
				if not name:
					name = codes.get((fix_bc1,fix_bc2))
					if not name:
						#finally give up
						unsorted += 1
						outfiles['unsorted'].write("".join(rec))
						continue
		outfiles[name].write("".join(rec))
	f.close()

print "%i total %.3f unsorted" % (total, unsorted/float(total))

for name in outfiles:
	outfiles[name].close()





