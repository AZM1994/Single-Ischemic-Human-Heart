#!/usr/bin/env python

import re
import sys
import random
import string

def main():

	CHROMLIST = open("chromlist.txt", "r")

	chromlist = []
	for x in CHROMLIST:
		if x.rstrip().find("_") == -1:
			chromlist.append(x.rstrip())

	INFILE = False
	OUTFILE = open("/filepath/sequence.part.0.k50.txt", "w")
	outfilecount = 0
	counter = 0
		
	for chrFile in chromlist:
		a = chrFile.split("/")
		b = a[len(a) - 1]
		thisChrom = b.split(".")[0]
		
		if INFILE:
			INFILE.close()
		INFILE = open(chrFile, "r")
		chr = []
		x = ""
		y = ""
		INFILE.readline()
		for x in INFILE:
			chr.append(x.rstrip())
		x = "".join(chr)
		y = x.upper()
		print "after read " + thisChrom
		sys.stdout.flush()
	
		for i in range(len(y) - 50):
			counter += 1
			thisRead = y[i:(i+50)]
			
			if counter % 150000000 == 0:   #### 150 million reads
				OUTFILE.close()	
				outfilecount += 1
				OUTFILE = open("/filepath/sequence.part." + str(outfilecount) + ".k50.txt", "w")
			
			OUTFILE.write("@" + thisChrom + "." + str(i))
			OUTFILE.write("\n")
			OUTFILE.write(thisRead)
			OUTFILE.write("\n")
			OUTFILE.write("+" + thisChrom + "." + str(i))
			OUTFILE.write("\n")
			OUTFILE.write("b" * len(thisRead))		
			OUTFILE.write("\n")

		INFILE.close()
	
	OUTFILE.close()


if __name__ == "__main__":
	main()
