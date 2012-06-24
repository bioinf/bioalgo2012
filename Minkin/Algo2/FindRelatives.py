import os
import re
import sys
import vcf

def parse23AndMe(path):    
    genotype = dict()  
    snpFile = open(path)
    for line in snpFile:        
        line = line.strip()
        if len(line) > 0 and line[0] != '#' and '-' not in line:
            line = line.split('\t')            
            genotype[line[0]] = line[3]                             
    snpFile.close()
    return genotype

def parseIllumina(path):
    genotype = dict()  
    snpFile = open(path)
    for line in snpFile:        
        line = line.strip()
        if len(line) > 0 and line[0] == '"' and '-' not in line:
            line = line.translate(None, '"').split(',')            
            genotype[line[0]] = line[3] 
    snpFile.close()    
    return genotype

def parseVCF(path):
    genotype = dict()
    snpFile = open(path)
    vcfReader = vcf.Reader(snpFile)    
    for record in vcfReader:
        gt = record.samples[0]['GT']   
        if gt == '0/0':
            alleles = record.REF + record.REF
        elif gt == '1/1':
            alleles = record.ALT[0] + record.ALT[0]
        else:
            alleles = record.ALT[0] + record.REF
        genotype[record.ID] = alleles   
    snpFile.close()
    return genotype

def parseFile(dir, entry):
    path = os.path.join(dir, entry)
    parser = [('vcf', parseVCF), ('23andme', parse23AndMe), ('illumina', parseIllumina)]
    for p in parser:
        if re.search(p[0], entry) != None:
            user = re.findall('user\d+', entry)
            if len(user) == 1:      
                user = re.findall('\d+', user[0])[0]
                return (user, p[1](path))
            return None    
    return None

def compareGenomes(genome1, prob, genome2):
    count = 0.0        
    for key in genome1:        
        if key in genome2 and genome1[key] == genome2[key]:
            p = 0.5    
            if key in prob:
                p = prob[key]
            count += p           
    return count

def SNPFrequency(record, query):
    allele = ''.join(record.alleles)    
    if len(allele) != 2:
        return None  
    for ch in query:
        if ch not in allele:
            return None        
    freq = 'AF'
    if freq not in record.INFO:
            return None    
    p = 1.0
    refCount = 0        
    for ch in query:
        pNow = float(record.INFO[freq])
        if allele.index(ch) != 1:
            pNow = 1 - pNow
            refCount += 1
        p *= pNow
    if len(query) == 2 and refCount == 1:
        p *= 2.0    
    return p

def calcProbability(genome, vcfFile):
    ret = dict()
    handle = open(vcfFile)
    vcfReader = vcf.Reader(handle)
    count = 1    
    MOD = 10000
    for record in vcfReader:      
        if count >= len(genome):
            break
        if record.ID in genome:
            count += 1
            if count % MOD == 0:
                print >> sys.stderr, 'Loaded: ', count
            prob = SNPFrequency(record, genome[record.ID])
            if prob != None: 
                ret[record.ID] = prob
    handle.close()
    return ret

if len(sys.argv) != 4:
    print "Usage: <genome file> <dir with other genomes> <vcf file>"
    exit(0)

count = 0
relative = []
fileDir = sys.argv[2]
genomeFile = sys.argv[1]
myGenome = parse23AndMe(genomeFile)
prob = calcProbability(myGenome, sys.argv[3])

for entry in os.listdir(fileDir):
    count += 1
    print >> sys.stderr, count, ' ', entry
    record = parseFile(fileDir, entry)        
    if record != None:
        score = compareGenomes(myGenome, prob, record[1])
        relative.append((score, record[0]))      
    
relative.sort(reverse=True)
print 'Closest relatives:'
print 'User\tScore'
for i in range(0, min(10, len(relative))):
    print relative[i][1], '\t', relative[i][0]