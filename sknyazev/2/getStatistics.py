#! /usr/bin/env python

import os
import sys
from optparse import OptionParser
import vcf
import _23andMe
from math import *
import numpy as np

class StatModes:
	ALT, REF = range(2)
class StatLogFields:
	AF, ASN_AF, AMR_AF, AFR_AF, EUR_AF = range(5)
class StatDataFields:
	ASN_AF, AMR_AF, AFR_AF, EUR_AF = range(4)
class Statistics:
	dataSet = []
	answers = []
	logs = [0,0,0,0,0]
	valueThreshold = 0.01
	deltaThreshold = 0.0
	smooth = 10000
	@staticmethod
	def smoothData(data):
		oldSum = 0
		newSum = 0
		for i in range(len(data)):
			oldSum += data[i]
			data[i] += Statistics.smooth
			newSum += data[i]
		for i in range(len(data)):
			data[i] = float(data[i]) * (float(oldSum)/float(newSum))
		return data
	@staticmethod
	def goodSnp(info, mode):
		values = []
		infoFields = ["AF", "ASN_AF", "AMR_AF", "AFR_AF", "EUR_AF"]
		for key in infoFields:
			if key not in info:
				return False
			if mode == "ALT":
				value = info[key]
			else:
				value = 1.0 - info[key]
			if  value < Statistics.valueThreshold:
				return False
			values.append(value)
		if max(values) - min(values) > Statistics.deltaThreshold:
			return True
		return False
	@staticmethod
	def getData(info, mode):
		dataFields = ["ASN_AF", "AMR_AF", "AFR_AF", "EUR_AF"]
		data = []
		if mode == "ALT":
			for key in dataFields:
				data.append(info[key])
#			answer = info["AF"]
		else:
			for key in dataFields:
				data.append(1 - info[key])
#			answer = 1 - info["AF"]
		data = Statistics.smoothData(data)
		answer = float(sum(data))/float(len(data))
		return data, answer

	def update(self, info, mode):
		data, answer = self.getData(info, mode)
		self.dataSet.append(data)
		self.answers.append(answer)
		self.updateLogs(data,answer)
	def updateLogs(self, data, answer):
		self.logs[StatLogFields.AF] += log(answer)
		self.logs[StatLogFields.ASN_AF] += log(data[StatDataFields.ASN_AF])
		self.logs[StatLogFields.AMR_AF] += log(data[StatDataFields.AMR_AF])
		self.logs[StatLogFields.AFR_AF] += log(data[StatDataFields.AFR_AF])
		self.logs[StatLogFields.EUR_AF] += log(data[StatDataFields.EUR_AF])

usage = "usage: %prog [options]"
description = "Description: Get snps race statistics."
parser = OptionParser(usage=usage, description=description)
parser.add_option("-f","--vcf", dest="vcfPath", help="vcf database FILE", metavar = "FILE")
parser.add_option("-g","--_23andme", dest="_23andMePath", help="23andMe database FILE", metavar = "FILE")
parser.add_option("-o", "--output_file", dest="outputPath", help="statistics FILE", metavar ="FILE")
(options, args) = parser.parse_args(sys.argv)

if options.vcfPath is None or options._23andMePath is None:
	parser.print_help()
	exit()

if options.outputPath is None:
	path = os.path.splitext(os.path.basename(options._23andMePath))
	options.outputPath = path[0] + "_statistics.txt"

snps = _23andMe._23andMe.loadSnps(options._23andMePath)
vcf_reader = vcf.VCFReader(open(options.vcfPath, 'rb'))

statistics = Statistics()

for record in vcf_reader:
	if record.ID in snps:
#		print snps[record.ID]
		if len(snps[record.ID]) == 2:
			if snps[record.ID][0] == '-' or snps[record.ID][1] == '-':
				continue
			if(snps[record.ID][0] == record.REF and snps[record.ID][1] == record.REF):
				if Statistics.goodSnp(record.INFO, "REF"):
					statistics.update(record.INFO, "REF")
					statistics.update(record.INFO, "REF")
			if((snps[record.ID][0] == record.REF and snps[record.ID][1] != record.REF) or
			(snps[record.ID][0] != record.REF and snps[record.ID][1] == record.REF)):
				if (Statistics.goodSnp(record.INFO, "REF")
				and Statistics.goodSnp(record.INFO, "ALT")):
					statistics.update(record.INFO, "REF")
					statistics.update(record.INFO, "ALT")
					statistics.updateLogs([2,2,2,2],2)
			if(snps[record.ID][0] != record.REF and snps[record.ID][1] != record.REF):
				if Statistics.goodSnp(record.INFO, "ALT"):
					statistics.update(record.INFO, "ALT")
					statistics.update(record.INFO, "ALT")
		elif len(snps[record.ID]) == 1:
			if(snps[record.ID] == record.REF):
				if Statistics.goodSnp(record.INFO, "REF"):
					statistics.update(record.INFO, "REF")
			else:
				if Statistics.goodSnp(record.INFO, "ALT"):
					statistics.update(record.INFO, "ALT")
		else:
			print snps[record.ID]

raceProbs = np.linalg.lstsq(statistics.dataSet, statistics.answers)
print raceProbs[0]

commonProbLog = 0
for i in range(4):
	probLog = (statistics.logs[i+1] + log(raceProbs[0][i]))
	if i == 0:
		commonProbLog = probLog
		continue
	if probLog - commonProbLog > 20:
		y = probLog - commonProbLog
	else:
		y = log(1 + exp(probLog - commonProbLog))
	commonProbLog = commonProbLog + y
results = []
results.append(str(exp(statistics.logs[StatLogFields.ASN_AF] + log(raceProbs[0][StatDataFields.ASN_AF])
- commonProbLog)))
results.append(str(exp(statistics.logs[StatLogFields.AMR_AF] + log(raceProbs[0][StatDataFields.AMR_AF])
- commonProbLog)))
results.append(str(exp(statistics.logs[StatLogFields.AFR_AF] + log(raceProbs[0][StatDataFields.AFR_AF])
- commonProbLog)))
results.append(str(exp(statistics.logs[StatLogFields.EUR_AF] + log(raceProbs[0][StatDataFields.EUR_AF])
- commonProbLog)))
f = open(options.outputPath, 'w')
f.write("#Probabilites:\n")
f.write("#Asian\t\tAmerican\tAfrican\t\tEuropean\n")

for i in range(4):
	f.write(str(results[i]))
	if i != 3:
		f.write('\t')
	else:
		f.write('\n')
f.close()
