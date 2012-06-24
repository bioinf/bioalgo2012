import sys
import vcf
from numpy import *

#VCF lib can be found here: https://github.com/jdoughertyii/PyVCF

globalF = 'AF'
asianF = 'ASN_AF'
americanF = 'AMR_AF'
africanF = 'AFR_AF'
europeanF = 'EUR_AF'
freqName = ['asian    ',
            'american ',
            'african  ',
            'european ',
            'global   ']
freqList = [asianF, americanF, africanF, europeanF, globalF]
 
def calcProbability(record, query):
    allele = ''.join(record.alleles)    
    if len(allele) != 2:
        return None  
    for ch in query:
        if ch not in allele:
            return None        
    ret = []
    for freq in freqList:
        p = 1.0
        refCount = 0
        if freq not in record.INFO:
            return None
        for ch in query:
            pNow = float(record.INFO[freq])
            if allele.index(ch) != 1:
                pNow = 1 - pNow
                refCount += 1
            p *= pNow
        if len(query) == 2 and refCount == 1:
            p *= 2.0    
        if p < 1E-10:
            return None    
        ret.append(p)     
    return ret

def getEstimationData(record):     
    ret = []
    for freqValue in freqList:
        if freqValue not in record.INFO:
            return None
        ret.append(record.INFO[freqValue])
    return ret

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print "Usage: <vcf> <23andme file> [cutoff]"
    exit(0)

cutoff = -1
if len(sys.argv) == 4:
    cutoff = int(sys.argv[3])

MOD = 100000
genotype = dict()
snpFile = open(sys.argv[2])
for line in snpFile:
    if len(genotype) % MOD == 0 and len(genotype) > 0:
        print >> sys.stderr, 'Loaded: ', len(genotype)
    line = line.strip()
    if len(line) > 0 and line[0] != '#':
        line = line.split('\t')        
        genotype[line[0]] = line[3]  

a = []
b = []
logFreq = [0.0] * len(freqList)
total = 0
counter = 0
vcfReader = vcf.Reader(open(sys.argv[1]))
for record in vcfReader:
    counter += 1
    if counter % MOD == 0:
        print >> sys.stderr, "Processed: ", counter    
    estData = getEstimationData(record)
    if estData != None:        
        a.append(estData[0:4])
        b.append(estData[4]) 
    if record.ID in genotype:
        ret = calcProbability(record, genotype[record.ID][0:1])
        if ret != None:
            total += 1              
            logFreq = [logFreq[i] + log(ret[i]) for i in range(0, len(ret))]                                      
            if total >= cutoff and cutoff > 0:
                break

print 'Total SNP considered = ', total
solution = linalg.lstsq(a, b)
print 'Estimated hypothesis:'
for i in range(0, 4):
    print freqName[i], ' -- ', solution[0][i]

acc = logFreq[0] + log(solution[0][0])
for i in range(1, len(logFreq) - 1):
    acc = logaddexp(acc, logFreq[i] + log(solution[0][i]))

print 'Probability that ' + sys.argv[2] + ' belongs to:'    
for i in range(0, 4):
    print freqName[i], ' -- ', exp(logFreq[i] - acc) * solution[0][i]