import re

class _23andMe:
	class Fields:
		ID, CHR, POS, SNP = range(4)
	@staticmethod
	def loadSnps(path):
		snps = {}
		for line in open(path, 'r'):
			if re.match('#', line):
				continue
			record = re.split('[\t\n\r]',line)
			if len(record) < 4:
				continue
			snps[record[_23andMe.Fields.ID]] = record[_23andMe.Fields.SNP]
		return snps
