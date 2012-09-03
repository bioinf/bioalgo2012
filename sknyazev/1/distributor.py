#! /usr/bin/python
import os

names = [ str(i) for i in range ( 1, 23 ) ] 
names.append("X")
names.append("Y")

###
### cuts regions from chromosomes that correspond to sequences
###  
def sequence( path ) :
                  
    content = ""
    if ( os.path.exists( path ) ) :
        for line in open ( path ) :
            if ">" in line :
                continue
            else :
                content += line.strip() 
    return content
 
### 
### distributes regions into files according to 
### their class membership   
###

def distribute( pathcls, pathchr, pathout ) :
    
        classes = []
        files = []
        positions = []
        for i in names : 
    
            infile_cls = pathcls + "chr" + i
            infile_chr = pathchr + "chr" + i + ".fasta" 
            if ( os.path.exists(infile_cls) and os.path.exists( infile_chr ) ) :
                print i
                content = sequence ( infile_chr ) 
                for line in open( infile_cls ) :
                    
                    res = line.split(" ")
		    if len(res) < 3:
			res.append("0")
                    if not (  int( res[2] ) in classes ) :
                        files.append( open ( (pathout + res[2]).strip(), 'w+' ) )
                        classes.append( int( res[2] ) )
                        positions.append( 0 )
                    
                    ind = classes.index( int( res[2] ) )
                    positions[ ind ] += 1
                    files[ ind ].write( ">" + str( positions[ ind ] ) \
		    + " chr" + str(i) + " " + str(int(res[0])) + " " + str(int(res[1])) + "\n" )
                    files[ ind ].write( content[ int( res[0] ) : int( res[1]) + 1 ] + "\n" )

        for k in files :
            k.close()




# example :
# distribute("/home/ksenia/Documents/bioalgo/train/","/home/ksenia/Documents/bioalgo/chrs/", "/home/ksenia/Documents/bioalgo/classes/" )
