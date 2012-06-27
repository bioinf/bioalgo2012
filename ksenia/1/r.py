#! /usr/bin/python

from aligner import gc, seq_range, header, viermers
from distributor import sequence
import operator

def kmers_input( seq, patterns )   :
	common = seq_range( seq, patterns, 20 )
	d = [] 
	for i in common :
		d.append(i[0])
	#debug print d

	p = patterns.keys()
	res = []
	#debug print "common: ", common
	#debug print "p: ", p
	#debug exit()
	for i in range( len(p) ) :
		if ( p[i] in d ) :
			res.append(1)
		else :
			res.append(0)
	return res

def output_classes( path, limit ) :

	f = open(path,"w+")

	for i in range(1,138) :
		for j in range(limit) :
			f.write( str(i) + "," )
			
	f.close()
			

def kmers_write( path, ar ) :
	for i in range( len(ar) ) :
		f = open( path+str(i), "a+" )
		#debug print "id: ", i, ar[i]
		f.write(str( ar[i] ) + ",")
		f.close()
		
def test( path_classes, path_out, limit ) :
	chrs = [str(i) for i in range(1,23) ]
        chrs.append('X')
        chrs.append('Y')
	v = viermers()
	

	for i in chrs :
			
		f = open( path_classes + str(i) )
		l = limit
		for seq in f :	
			if ( l == 0 ) : 
				break
			if not header( seq ) :
				res = kmers_input( seq, v )
				kmers_write( path_out, res )
				
				l -= 1
		f.close()		

def train( path_classes, path_out, limit, output_class=True ) :
	v = viermers()
	
	cl = open("/labnas/students/ksenia/bioalgo/r/kmers_test/classes","w")

	for i in range( 1,138 ) :
			
		f = open( path_classes + str(i) )
		l = limit
		for seq in f :	
			if ( l == 0 ) : 
				break
			if not header( seq ) :
				if ( output_class ) :
					cl.write( str(i) + "," )
				res = kmers_input( seq, v )
				kmers_write( path_out, res )
				#print "l: ", l
				
				l -= 1
		#if ( l != 0 ) :
		#	print "too few sequences: ", l
		#	exit (1)
		f.close()		
	cl.close()

def seq_features( seq, classes, gcf, length ) :
	classes.write( str(i) + "," )
	gcf.write( str( gc(seq) ) + "," )
	length.write( str( len(seq) ) + "," ) 

def features( file_names, path_classes, path_out ) :
	classes = open( path_out + "classes", "w+" )
	gcf = open( path_out + "gc", "w+" )
	length = open( path_out + "length", "w+" )
	#kmers = open( path_out + "kmers", "w+" )

	patterns = viermers()
	
	for i in file_names :
		print "class ", i
		f = open(path_classes + i, "r")
		lim = 20	
		for seq in f :
			if not header( seq ) :
				seq_features( seq, classes, gcf, length )
				#common = seq_range( seq, patterns, 10 ) 
				#for el in common :
				#    kmers.write(el + " ")
				#kmers.write("\n")
				lim -= 1
		if lim < 0 :
		    break
		f.close()
	classes.close()
	gcf.close()
	length.close()
	#kmers.close()

def convert( path ) :
# path - file of model prediction

	chrs = [str(i) for i in range(1,23) ]
	chrs.append('X')
	chrs.append('Y')

	res = open( path, "r" )

	counter = 0
	path_out = "./tested_final/"
	for i in chrs :
		print "chr" + i
		f_out = open( path_out + "chr" + i , "w" )
		for line in open( "/labnas/bioalgo/test/chr" + i, "r" ) :	

			counter += 1
			#print line
			r = res.readline()
			r = r.strip()
			if r == "" :
				r = "1"
			line = line.strip()
		#	print line + " " + r 
			f_out.write( line + " " + r + "\n")

		f_out.close()

	print counter

	res.close()

	
if __name__ == "__main__" :

	convert("./predicted")

