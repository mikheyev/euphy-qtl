# Input: fastq reads demultiplexed
# Output: fastq reads with 4-bases prefix added. Ready to remove duplicates.

import sys,gzip,itertools,glob,pdb,os

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

pathFile=str(sys.argv[1])
##Obtain basename and dirname for output file
fileName=os.path.basename(pathFile)
fileName=fileName.split('.')[0]
dirName=os.path.dirname(pathFile)

#out = open(dirName+"/prefix/"+fileName+".fq", "w")
out=open(str(sys.argv[2]),"w")

f = open(pathFile)
for rec in parse_fastq(f):
	barC=rec[0].split(':')[9]
	barC = barC.replace("\n","")
	#rec[1]=barC[7:-7]+rec[1]
	#batch3
	#rec[1]=barC[7:-9]+rec[1]
	#batch4-batch8
	rec[1]=barC[7:-10]+rec[1]
	rec[3]="IIII"+rec[3]
	barC=rec[0].split(':')[9]
	out.write("".join(rec))
out.close()

