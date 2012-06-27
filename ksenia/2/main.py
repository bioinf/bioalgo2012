from pysam import VCF
from math import log
from math import exp
from string import find

import re
import csv
import os

def loadGlobalData () :
	snps = {}
	haplotype = {}

	counter = 0
	for line in open("/labnas/students/slebedev/opensnp/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf", "r") :
		if ( not ("##") in line ) and ( not ("#CHROM") in line ) :
			if counter > 700000 :
				return snps, haplotype
			res = re.split( r'[\t,]', line )
			if len( res ) == 8 :
			
				info = res[-1].split(";")
				l = []
				for i in info :
					if "AF" in i :
						if not "LDAF" in i:
							
							l.append( i.strip() )		
				#print l
				snps[ res[1] ] = l 
				haplotype[ res[1] ] = res[4]
				#haplotype.append( res[4] )
				#print "(", res[1] , l, ")", res[4]
				counter += 1
	return snps, haplotype
	

def loadSNP ( path ) :
	
	snps = []
	haplotype = []
	counter = 0
	for line in open ( path ) :
		#if counter > 500 :
		#	break
		if not "#" in line :
			
			res = re.split( r'[\t\n,]', line ) #re.findall( r'\w+', line )
			if len(res) > 3 :
				counter += 1
				snps.append( res[2].strip() )
				haplotype.append( res[3].strip() )
			#else :
				#print "size of data is too small!"
				#print res
	
	return snps, haplotype
			
def search( l, key ) :

	#print "---------"	
	i = 0
	for a,b in l:
		#print a,b
		if str(a) == str(key) :
			print "match"
			return b, i
		i += 1
	return None, None

def count ( own_snps, ho, global_snps, hg ) :

	print "counting probabilities..."	

	asn = 0; afr = 0; amr = 0; eur = 0; 
	epsilon = 0.0001	
	 
	i = -1
	for snp in own_snps :
		i += 1
		if snp in global_snps : 
			r = global_snps[ snp ]
			h = ho[i]
			f = find( h,hg[snp] )

			if r != None and f != -1 : 

			#	print snp, h, hg[snp]
				names = []
				vals = []
				for el in r :
					pair = el.split("=")
					names.append( pair[0] )
					vals.append( pair[1] )
				
				if "ASN_AF" in names :
					asn += log (  1 - float( vals[names.index("ASN_AF")]) + epsilon )  
				else :
					asn += log (  epsilon )  
				if "EUR_AF" in names :
					eur += log (  1 - float( vals[names.index("EUR_AF")]) + epsilon )  
				else :
					eur += log (  epsilon )  
				if "AMR_AF" in names :
					amr += log (  1 - float( vals[names.index("AMR_AF")]) + epsilon )  
				else :
					amr += log (  epsilon )  
				if "AFR_AF" in names :
					afr += log ( 1 - float( vals[names.index("AFR_AF")]) + epsilon )  
				else :
					afr += log (  epsilon )  
	
	print asn, eur, afr, amr
	m = max ( [ asn, eur, afr, amr ] )
	asn = exp( asn - m ); eur = exp ( eur - m ); afr = exp ( afr - m ); amr = exp ( amr - m )
	print asn, eur, afr, amr
	s =  sum ( [ asn, afr, amr, eur ] ) 
	print "s:",
	p_afr = ( afr / s )
	p_asn = ( asn / s )
	p_eur = ( eur / s)
	p_amr = ( amr / s )

	return [ p_asn, p_eur, p_afr, p_amr ]	

def initDict ()  :

	dictSNP = ()
	for i in [ "A", "C", "G", "T" ] :
		for j in [ "A", "C", "G", "T" ] :

			dictSNP [ i+j ] = k

	return dictSNP

def similarity ( snps1, h1, snps2, h2, global_snps ) :

	#print "inside of similarity"	
	s1 = set( zip( snps1, h1 ) )
	s2 = set( zip( snps2, h2 ) )
	
	common = s1.intersection(s2)
	score = 0
	for el in common :
		r = global_snps.get ( el[0] )
#		print "el: ", el[0], el[1] 
#		print "r: ", r
		if r :
			for p in r :
				pair = p.split("=")
				if ( p[0] == "AF" ) :
					score += 1 - float( p[1] )

	return score

def find_max( ar, pos ) :
	cut = [ el[pos] for el in ar ]
	return max( ( set( cut )).difference( set(["-"]) ) , key=cut.count)


def phenotype( sims, path ) :
	
	ns = []
	for item in sims :
		r = item[1].split("_")	
		ns.append( re.findall( r'\d+', r[0]) )

	r = csv.reader(open(path,"rb"), delimiter=";")
	ph = []
	for row in r :
		id = re.findall(r"^\d+$", row[0] )
		if id and id in ns :
			ph.append( row  )	
	res = []
	for pos in range(3,14) :
		res.append( find_max( ph,pos ) )
	return res

### task 1
print "loading global data..."
global_snps, global_haplotype = loadGlobalData()
print "global snps: ", len(global_snps)

print "loading local snps..."
local_snps, haplotype = loadSNP( "/labnas/students/ksenia/algo/2/genome_Ksenia_Krasheninnikova_Full_20120530043820.txt" )
print "local snps: ", len(local_snps)

print "processing..."
print count ( local_snps, haplotype, global_snps, global_haplotype )
#print " order of probabilities: african, european, asian, american:"
### end task 1

### task2
##path = "/labnas/students/slebedev/opensnp/opensnp/"
##files = os.listdir( path ) 
##print "task2"

##sims = []
##i = 0

##for file in files :

#debug	print file
##	if ( file.split(".")[-1] == "txt" ) :
##		neighbour_snps, neighbour_h = loadSNP( path+file ) 
##		sims.append( ( similarity( local_snps, haplotype, neighbour_snps, neighbour_h , global_snps ), files[i] ) ) 
##	i += 1
#debug	if i > 11 :
#debug		break

##sims.sort( key=( lambda a : a[0] )  )
#print sims
#print "it was sims"
#debug file = open("./result_2_2.txt" , "w+")
#debug for item in sims[:10] :
#debug  print>>file, item
#debug file.close()

#debug 
##print sims[:10]
### end task2

### task3
##print "task 3"
##print phenotype ( sims[:10], "/labnas/students/slebedev/opensnp/opensnp/phenotypes_201205080944.csv" )
### end task3

			
			
		
	
