import vcf
import math
import sys

##reader = vcf.Reader(file("data/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf"))


class Position:
    
    def __init__(self, reader, rec_to_pos_chrom):
        self._reader = reader
        self.chrom = "a"
        self.chrom_ind = -1
        self.converter = rec_to_pos_chrom

    def __lt__(self, other):
        return (self.chrom_ind < other.chrom_ind) or (self.chrom_ind == other.chrom_ind and self.pos < other.pos)

    def __eq__(self, other):
        return self.chrom_ind == other.chrom_ind and self.pos == other.pos

    def __ne__(self, other):
        return not self.__eq__(other)


    def process_record(self, record):
        self.record = record
        chrom, pos = self.converter(self.record)  
        self.pos = int(pos)
        ##TODO: change calculating of chrom_ind
        if self.chrom <> chrom:
            self.chrom_ind += 1
            self.chrom = chrom

    
    def next(self):
        try:
            self.process_record(self._reader.next())
            return True
        except StopIteration:
            return False



def gen_converter(line):
    return line.split('\t')[1:3]

def vcf_converter(vcf_record):
    return vcf_record.CHROM, vcf_record.POS

vcf_reader = vcf.Reader(file(sys.argv[1]))
reader_23andme = open(sys.argv[2])
line = reader_23andme.next()
while line.startswith('#'):
    line = reader_23andme.next()



vcf_position = Position(vcf_reader, vcf_converter)
vcf_position.next()

gen_position = Position(reader_23andme, gen_converter) 
gen_position.process_record(line)


eps = 0.00001
amr_pr = math.log(0.25)
afr_pr = math.log(0.25)
eur_pr = math.log(0.25)
asn_pr = math.log(0.25)

eof_reached = False


def alt_prob_calc(prob, race_af, af):
    return prob + math.log( race_af / af)

def ref_prob_calc(prob, race_af, af):
    return  prob + math.log( (1.0 - race_af) / (1.0 - af))




while True:
    
    while (vcf_position <> gen_position):
       # print gen_chrom , " " , gen_pos , " " , vcf_chrom , " " , vcf_pos
        if (gen_position < vcf_position ):
          #  print "first"            
            if not gen_position.next():
                eof_reached = True
                break     
        else:
            # print "second"
            if not vcf_position.next():
                eof_reached = False
                break

    if eof_reached:
        break

    genome_gen = gen_position.record.split('\t')[3]

    vcf_record = vcf_position.record

   
    info = vcf_record.INFO
    
    af = info.get("AF")
    
        
    amr_af = info.get("AMR_AF", eps)
    afr_af = info.get("AFR_AF", eps)
    eur_af = info.get("EUR_AF", eps)
    asn_af = info.get("ASN_AF", eps)
   
    lst = [amr_af, afr_af, eur_af, asn_af]
#    undef_c = len([i for i in lst if i == -1])
#    if undef_c <> 0:
#        sum_def = sum([i for i in lst if i <> -1]) 
#        prob = (4 * af - sum_def ) / undef_c     
#        for i in range(4):
#            if lst[i] == -1:
#                lst[i] = prob

#    lst = [i * 0.25 for i in lst]
    if (af):
        #print vcf_record.REF, vcf_record.ALT , genome_gen  
        #print info        
        """    if vcf_record.REF == genome_gen[0]:
                print "ref found"
                if amr_af > 0:
                    amr_pr = ref_prob_calc(amr_pr, amr_af, af) 
                if afr_af > 0:
                    afr_pr = ref_prob_calc(afr_pr, afr_af, af)   
                
                if eur_af > 0:                
                    eur_pr = ref_prob_calc(eur_pr, eur_af, af)
                    
                if asn_af > 0:    
                    asn_pr = ref_prob_calc(asn_pr, afr_af, af) 

                if not gen_position.next():
                    break
          """
            
        for alt in vcf_record.ALT:
            
            if (str(alt)[0] == genome_gen[0]):
               #print "alt found"
                
                amr_pr = alt_prob_calc(amr_pr, lst[0], af)
            
                afr_pr = alt_prob_calc(afr_pr, lst[1], af)
            
                eur_pr = alt_prob_calc(eur_pr, lst[2], af)
            
                asn_pr = alt_prob_calc(asn_pr, lst[3], af)
                break

        print "AMR ", amr_pr#, math.exp(amr_pr)
        print "AFR ", afr_pr#, math.exp(afr_pr)
        print "EUR ", eur_pr#, math.exp(eur_pr)
        print "ASN ", asn_pr, "\n" #, math.exp(asn_pr)

    if not gen_position.next():
        break

   ##print reader.next().CHROM
#record = vcf_reader.next()
#while True:
#    count += 1
#    try:
#        record = reader.next()
#    except StopIteration:
#        break
    
#irint count
   
