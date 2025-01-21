#!/usr/bin/env python

import sys


def main():

	infilename = sys.argv[1]
	outfilename = sys.argv[2]
	statfilename = sys.argv[3]
	binfilename = sys.argv[4]
	chromsizefile = sys.argv[5] # "hg19.chrom.sizes.txt"
	#binfilename = "hg19.bin.boundaries.50k.bowtie.k50.sorted.txt"

	chrominfo = fileToDictionary(chromsizefile, 0)
	bins = fileToArray(binfilename, 0)
	INFILE = open(infilename, "r")
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	binCounts = []
	for i in range(len(bins)):
		binCounts.append(0)

	print len(binCounts)
	print len(bins)

	counter = 0
	totalReads = 0
	prevChrompos = ""
	for x in INFILE:
		if x[0] != "@":
			arow = x.rstrip().split("\t")
			#print arow
			thisChrom = "chr"+arow[2]
			thisChrompos = arow[3]
			#print thisChrom
			if thisChrom.find("_") > -1:
				#print thisChrom
				continue
			if thisChrom == "chrM":
				#print thisChrom
				continue
			if thisChrom == "":
				continue
			if chrominfo.has_key(thisChrom):
				pass
			else:
				continue

			totalReads += 1
				
			thisChrominfo = chrominfo[thisChrom]
			thisAbspos = long(thisChrompos) + long(thisChrominfo[2])
			
			counter += 1
			#print counter
			
			indexUp = len(bins) - 1
			indexDown = 0
			indexMid = int((indexUp - indexDown) / 2.0)

			while True:
				if thisAbspos >= long(bins[indexMid][2]):
					indexDown = indexMid + 0
					indexMid = int((indexUp - indexDown) / 2.0) + indexMid
				else:
					indexUp = indexMid + 0
					indexMid = int((indexUp - indexDown) / 2.0) + indexDown

				if indexUp - indexDown < 2:
					break

			binCounts[indexDown] += 1
			prevChrompos = thisChrompos
	
	print len(bins)
	print counter
	'''
	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / (float(counter) / float(len(bins)))
		OUTFILE.write("\t".join(bins[i][0:3]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(thisRatio))
		OUTFILE.write("\n")
	'''

	cont = ""
	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / (float(counter) / float(len(bins)))
		cont += "\t".join(bins[i][0:3])
		cont += "\t"
		cont += str(binCounts[i])
		cont += "\t"
		cont += str(thisRatio)
		cont += "\n"
	OUTFILE.write(cont)
	binCounts.sort()
	
	STATFILE.write("TotalReads\tMedianBinCount\n")
	STATFILE.write(str(totalReads))
	STATFILE.write("\t")
	STATFILE.write(str(binCounts[len(bins)/2]))
	STATFILE.write("\n")

	INFILE.close()
	OUTFILE.close()
	STATFILE.close()


def fileToDictionary(inputFile, indexColumn):
	input = open(inputFile, "r")

	rd = dict()
#	input.readline()
	for x in input:
		arow = x.rstrip().split("\t")
		id = arow[indexColumn]
		if rd.has_key(id):
			#rd[id].append(arow)
			print "duplicate knowngene id = " + id
			print "arow =   " + str(arow)
			print "rd[id] = " + str(rd[id])
		else:
			rd[id] = arow
		
	input.close()
	return(rd)


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []

	for i in range(skipFirst):
		input.readline()

	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)
		
	input.close()
	return(ra)


if __name__ == "__main__":
	main()
