#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import bisect


def main():

	#bincount = 50000
	MAP = open("hg19.chrom.mappable.bowtie.k50.txt", "r")
	GOOD = open("hg19.goodzones.bowtie.k50.bed", "r")
	CHROMLEN = open("hg19.chrom.sizes.txt", "r")
	#OUTFILE = open("hg19.bin.boundaries.50k.bowtie.k50.txt", "w")


	#bincount = 6000
	#OUTFILE = open("hg19.bin.boundaries.6k.bowtie.k50.txt", "w")	

	bincount = 20000
	OUTFILE = open("hg19.bin.boundaries.20k.bowtie.k50.txt", "w")

	
	#OUTFILE.write("chrom\tbin.start.chrompos\tbin.start.abspos\tbin.end.chrompos\tbin.length\tmappable.positions\n")

	chromlen = dict()
	for x in CHROMLEN:
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisChromlen = int(arow[1])
		thisAbspos = long(arow[2])
		chromlen[thisChrom] = [thisChromlen, thisAbspos]

	chromarray = []
	chroms = dict()
	totalLength = long(0)
	for x in MAP:
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisLength = long(arow[1])
		if thisChrom == "chrM":
			continue
		totalLength += thisLength
		chroms[thisChrom] = thisLength

	bincountUsed = 0
	for k, v in chroms.items():
		chromBincount = float(bincount) * (float(v) / float(totalLength))
		i = int(chromBincount)
		bincountUsed += i 
		r = chromBincount - i
		chromarray.append([k, i, r])

	a = []
	for i in chromarray:
		bisect.insort(a, (-i[2], i))

	chromarray = []
	for j in a:
		chromarray.append(j[1]) ##  This j[1] index has to match the index of i in the bisect.insort 2nd parameter.

	a = []
	
	print chromarray	

	remain = bincount - bincountUsed
	for i in range(remain):
		chromarray[i][1] += 1

	chroms2 = dict()
	for i in range(len(chromarray)):
		chromlength = chroms[chromarray[i][0]]
		chrombins = chromarray[i][1]
		binlength = float(chromlength) / float(chrombins)
		chroms2[chromarray[i][0]] = [chrombins, binlength]

	
	chromorder = sorted(chroms2.keys())
	print chromorder

	print 
	print "Starting to get bin boundaries"
	print
	goodEOF = False
	for chrom in chromorder:
		print chrom
		firstBin = True
		#  position GOOD file
		x = GOOD.readline()
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisStart = int(arow[1])
		thisEnd = int(arow[2])
		print "new chrom", arow
		while thisChrom != chrom:
			x = GOOD.readline()
			arow = x.rstrip().split("\t")
			if len(x) == 0:
				goodEOF = True
				break
			thisChrom = arow[0]
			thisStart = int(arow[1])
			thisEnd = int(arow[2])
			print "position chrom", arow
		if goodEOF:
			break
		print "after position"
		currentStart = thisStart 
		chromBincount = chroms2[chrom][0]
		chromBinlength = chroms2[chrom][1]
		chromExcess = chromBinlength - int(chromBinlength)
		currentExcess = 0.0
		thisBincount = 0
		binStart = 0
		binEnd = 0
		currentLength = 0
		print chromBincount, chromBinlength
		while thisBincount < chromBincount:
			thisBincount += 1
			print "thisBincount", thisBincount
			thisBinlength = int(chromBinlength)
			currentExcess += chromExcess
			print "currentExcess", currentExcess
			if currentExcess >= 1.0:
				currentExcess -= 1.0
				thisBinlength += 1
			print "thisBinlength", thisBinlength
			binStart = currentStart
			currentLength = 0
			if (binStart + thisBinlength) < thisEnd:
				print "got bin from current GOOD"
				binEnd = binStart + thisBinlength
				currentStart = binStart + thisBinlength
				currentLength = thisBinlength
			else:
				print "getting more GOOD"
				currentLength += (thisEnd - currentStart)
				print thisEnd, currentStart, currentLength, thisBinlength
				while currentLength < thisBinlength:
					x = GOOD.readline()
					print x
					if len(x) == 0:
						goodEOF = True
						break
					arow = x.rstrip().split("\t")
					thisChrom = arow[0]
					thisStart = int(arow[1])
					thisEnd = int(arow[2])
					print "adding length", arow
					if thisChrom != chrom:
						print "ERROR: Past end of chrom.", thisChrom, chrom
						break

					if (thisEnd - thisStart) < (thisBinlength - currentLength):
						currentLength += (thisEnd - thisStart)
					else:
						currentStart = thisStart + (thisBinlength - currentLength)
						currentLength = thisBinlength
					print "currentLength", currentLength
				binEnd = currentStart

			if firstBin:
				binStart = 0
			if thisBincount == chromBincount:
				binEnd = chromlen[chrom][0]

			binStartAbspos = chromlen[chrom][1] + binStart

			OUTFILE.write(chrom)
			OUTFILE.write("\t")
			OUTFILE.write(str(binStart))
			OUTFILE.write("\t")
			OUTFILE.write(str(binStartAbspos))
			OUTFILE.write("\t")
			OUTFILE.write(str(binEnd))
			OUTFILE.write("\t")
			OUTFILE.write(str(binEnd - binStart))
			OUTFILE.write("\t")
			OUTFILE.write(str(currentLength))
			OUTFILE.write("\n")

			firstBin = False

			if goodEOF:
				break


	OUTFILE.close()
	MAP.close()
	GOOD.close()


if __name__ == "__main__":
	main()
