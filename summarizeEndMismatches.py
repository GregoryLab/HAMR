#!/usr/local/bin/python
import sys

#version 1.0

########

alphabet = [',', '.', 'a', 'c', 'g', 't', 'A', 'C', 'G', 'T']
start_mismatch = dict((letter,int(0)) for letter in alphabet)
end_mismatch = dict((letter,int(0)) for letter in alphabet)
nonend_mismatch = dict((letter,int(0)) for letter in alphabet)

#######

inFile = open(sys.argv[1],'r')
lineNumber = 0
for line in inFile:
	lineNumber += 1
	lineArray=line.rstrip().split('\t')
	ntList=list(lineArray[4])
	count = 0
	while count < len(ntList):
		if ntList[count] == '$':
			count += 1
			end_mismatch[ntList[count]] += 1
		elif ntList[count] == '^':
			count += 2
			start_mismatch[ntList[count]] += 1
		elif ntList[count] in alphabet:
			nonend_mismatch[ntList[count]] += 1
		else:
			sys.exit("Unrecognized pileup symbol \""+str(ntList[count])+"\" at line "+str(lineNumber))
		count += 1
inFile.close()

#######

print "nt\tFirst\tMiddle\tLast\tPercentFirst\tPercentLast"
for letter in alphabet:
	total = float(start_mismatch[letter]+nonend_mismatch[letter]+end_mismatch[letter])
	print str(letter)+'\t'+str(start_mismatch[letter])+'\t'+str(nonend_mismatch[letter])+'\t'+str(end_mismatch[letter])+'\t'+str(start_mismatch[letter]/total*100)+'\t'+str(end_mismatch[letter]/total*100)

