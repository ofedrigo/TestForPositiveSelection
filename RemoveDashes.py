#! /usr/bin/env python
import optparse,sys,sets

def getFasta(fileName):
	fileHandle=open(fileName,"r")
	data=fileHandle.read().splitlines()
	fileHandle.close()
	results={}
	seq=""
	compName=""
	for line in data:
		if line[0]==">":
			if compName!="":
				results[compName]=seq
			compName=line[1:]
			seq=""
		else:
			seq=seq+line
	results[compName]=seq
	nsize=len(results[results.keys()[0]])
	ns=len(results.keys())
	return results,nsize,ns

def printFasta(fastaDict):
	for spec in fastaDict.keys():
		print ">"+spec
		print fastaDict[spec]

def removeIndels(fastaDict,nsize,ns):
	indels=sets.Set([])
	fastaDict2={}
	for spec in fastaDict.keys():
		for i in range(nsize):
			base=fastaDict[spec][i]
			if base.upper() not in ["A","C","G","T"]:
				indels.add(i)
	print indels
	nsize2=nsize-len(indels)
	for spec in fastaDict.keys():
		seq=""
		for i in range(nsize):
			if i not in list(indels):
				seq=seq+fastaDict[spec][i]
		fastaDict2[spec]=seq
	return fastaDict2,nsize2
	
def main():
	optParser = optparse.OptionParser(usage = "%prog input output",description="This script corrects Fasta alignments (remove indels and Ns)")
	 #epilog =  "Written by Olivier Fedrigo (ofedrigo@duke.edu)")
	if len( sys.argv ) <1:
		optParser.print_help()
		sys.exit(1)
	(opts, args) = optParser.parse_args()
	input=args[0]
	ouput=args[1]

	myData,nsize,ns=getFasta(input)
	myData2,nsize2=removeIndels(myData,nsize,ns)
	print "original size=",nsize,"bp"
	print "final size=",nsize2,"bp"
	
	outputHandle=open(ouput,"w")
	for i in range(len(myData2.keys())):
		spec=myData2.keys()[i]
		if i>0: outputHandle.write("\n")
		outputHandle.write(">"+spec+"\n")
		outputHandle.write(myData2[spec])
	outputHandle.close()
	
if __name__ == "__main__":
	main() 
