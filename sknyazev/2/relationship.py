#! /usr/bin/env python

import sys
import os
import _23andMe
import glob
import vcf
from optparse import OptionParser

def compareSnps(snp1,snp2):
	matches = 0
	mismatches = 0
	for i in snp1:
		if i == '-':
			return 0, 0
	for i in snp2:
		if i == '-':
			return 0, 0
	for i in snp1:
		if i in snp2:
			matches += 1
		else:
			mismatches += 1
	return matches, mismatches

def getCandidatesPaths(directory):
	currentDirectory = os.getcwd()
	os.chdir(directory)
	files = glob.glob('*.23andme.txt')
	os.chdir(currentDirectory)
	return files

def getWeight(snp, snpScores):
	if snp in snpScores:
		return snpScores[snp]
	else:
		return 0

usage = "usage: %prog [options]"
description = "Description: Find relationship."
parser = OptionParser(usage=usage, description=description)
parser.add_option("-f","--tested_person", dest="inputPath", help="23andMe FILE", metavar = "FILE")
parser.add_option("-v","--vcf", dest="vcfPath", help="vcf database FILE", metavar = "FILE")
parser.add_option("-d","--relationship_candidates_dir"
, dest="candidatesDir", help="candidates DIR", metavar = "DIR")
parser.add_option("-o", "--output_file", dest="outputPath", help="output FILE", metavar ="FILE")
(options, args) = parser.parse_args(sys.argv)

if options.inputPath is None or options.candidatesDir is None or options.vcfPath is None:
	parser.print_help()
	exit()

if options.candidatesDir[-1] != '/':
	options.candidatesDir = options.candidatesDir + "/"

if options.outputPath is None:
	path = os.path.splitext(os.path.basename(options.inputPath))
	options.outputPath = path[0] + "_relationship.txt"

files = getCandidatesPaths(options.candidatesDir)
if not files:
	print "No files in candidates directory"
	exit()
snps = _23andMe._23andMe.loadSnps(options.inputPath)
vcf_reader = vcf.VCFReader(open(options.vcfPath, 'rb'))
snpScores = {}
for record in vcf_reader:
	snpScores[record.ID] = record.INFO["AF"]

relations = []
for candidat in files:
	candidateSnps = _23andMe._23andMe.loadSnps(options.candidatesDir + candidat)
	snpsScore = 0
	snpsCount = 0
	for snp in candidateSnps:
		if snp in snps:
			matches, mismatches = compareSnps(snps[snp],candidateSnps[snp])
			snpsScore += matches * getWeight(snp, snpScores)
			snpsCount += (matches + mismatches)
	if snpsCount != 0:
		print (candidat + '\t' + str(float(snpsScore)/float(snpsCount)) + '\t' 
		+ str(snpsScore) + '\t' + str(snpsCount))
		relations.append((candidat, float(snpsScore)/float(snpsCount), snpsScore, snpsCount))
	else:
		relations.append((candidat, 0, 0, 0))
relations = sorted(relations, key = lambda x: x[1], reverse=True)
f = open(options.outputPath, 'w')
for record in relations:
	f.write(record[0] + '\t' + str(record[1]) + '\t' + str(record[2]) + '\t' + str(record[3]) + '\n')
f.close()

## usage ./relationship.py -d "/storage/labnas/students/slebedev/opensnp/raw/" -f genome_Sergey_Knyazev_Full_20120613051503.txt
