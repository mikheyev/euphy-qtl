import sys

files={}
with open(sys.argv[1]) as f:
	f.readline()
	for line in f:
		lineList=line.split("\t")
		if lineList[1] in files:
			files[lineList[1]]+=lineList[2]+"\t"+lineList[5][:7]+"\t"+lineList[5][-7:]+"\n"
		else:
			files[lineList[1]]=lineList[2]+"\t"+lineList[5][:7]+"\t"+lineList[5][-7:]+"\n"
	
conta=1		
for var in files:
	f = open('csvs/lane'+str(conta)+'_bcodes.csv', 'w')
	f.write(files[str(conta)])
	conta+=1
	f.close()

