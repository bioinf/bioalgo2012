#! /usr/bin/env python

import re
import sys
import os
from optparse import OptionParser
import _23andMe

class VCF:
	class Fields:
		CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = range(8)
	infoFields = ["AF", "ASN_AF", "AMR_AF", "AFR_AF", "EUR_AF"]
	valueThreshold = 0.1
	deltaThreshold = 0.1
	@staticmethod
	def loadSnps(path):
		snps = {}
		for line in open(path, 'r'):
			if re.match('#', line):
				continue
			record = re.split('[\t\n]',line)
			snps[record[VCF.Fields.ID]] = record[VCF.Fields.REF]
		return snps
	@staticmethod
	def getValue(info,key):
		pattern = re.compile(key + "=(?P<val>[0-9]*\.[0-9]+)")
		value = pattern.search(info)
		if value != None:
			return float(value.group('val'))
		return 0.0
	

usage = "usage: %prog [options]"
description = "Description: Filter vcf file by snps present at 23andMe file."
parser = OptionParser(usage=usage, description=description)
parser.add_option("-f","--vcf", dest="vcfPath", help="vcf FILE to be filtered", metavar = "FILE")
parser.add_option("-g","--_23andme", dest="_23andMePath", help="23andMe snps FILE", metavar = "FILE")
parser.add_option("-o", "--vcf_output", dest="vcfOutputPath", help="vcf output FILE", metavar ="FILE")
parser.add_option("-p", "--_23andme_output", dest="_23andMeOutputPath", help="23andMe output FILE"
, metavar ="FILE")
(options, args) = parser.parse_args(sys.argv)

if options.vcfPath is None or options._23andMePath is None:
	parser.print_help()
	exit()
if options.vcfOutputPath is None:
	path = os.path.splitext(os.path.basename(options.vcfPath))
	options.vcfOutputPath = path[0] + "_filtered" + path[1]
if options._23andMeOutputPath is None:
	path = os.path.splitext(os.path.basename(options._23andMePath))
	options._23andMeOutputPath = path[0] + "_filtered" + path[1]

# filter vcf
snps = _23andMe._23andMe.loadSnps(options._23andMePath)
f = open(options.vcfOutputPath, 'w')
for line in open(options.vcfPath, 'r'):
	if re.match("##", line) or re.match("#CHROM", line):
		f.write(line)
		continue
	record = re.split("[\t\n]", line)
	if record[VCF.Fields.ID] in snps:
			f.write(line)
f.close()

#filter snps
snps = VCF.loadSnps(options.vcfOutputPath)
f = open(options._23andMeOutputPath, 'w')
for line in open(options._23andMePath):
	if re.match("#", line):
		f.write(line)
		continue
	record = re.split("[\t\n]", line)
	if record[_23andMe._23andMe.Fields.ID] in snps:
		f.write(line)
