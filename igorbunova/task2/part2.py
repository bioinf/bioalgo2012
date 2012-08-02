import sys
import os

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

    def __get_chrom_ind(self, chrom_name):
        if (chrom_name == "X"):
            return 30
        elif (chrom_name == "Y"): 
            return 31
        elif (chrom_name == "MT"):
            return 32
        else:
            return int(chrom_name)
                
    def process_record(self, record):
        try:
            self.record = record
            chrom, pos = self.converter(self.record)  
            self.pos = int(pos)
            ##TODO: change calculating of chrom_ind
            if self.chrom <> chrom:
                self.chrom_ind = self.__get_chrom_ind(chrom)
                self.chrom = chrom
            return True
        except Exception:
            return False
            

    def next(self):
        try:
            while not self.process_record(self._reader.next()):
                pass
            return True
        except StopIteration:
            return False



def gen_converter(line):
    return line.split('\t')[1:3]


def skip_comments(reader):
    line = reader.next()
    while line.startswith('#'):
        line = reader.next()
    return line


def init_position(file_name):
    reader = open(file_name)
    line = skip_comments(reader)
    position = Position(reader, gen_converter)
    position.process_record(line)
   
    return position

def calc_inv_dist(file1, file2, result_file):
    
    pos1 = init_position(file1)
    pos2 = init_position(file2)
    eof_reached = False
    count = 0
    while True:
        
        while (pos1 <> pos2):
            increase_pos =  pos1 if pos1 < pos2 else pos2
            if (not increase_pos.next()):
                eof_reached = True
                break

        if (eof_reached):
            break

        gen1 = pos1.record.split('\t')[3]
        gen2 = pos2.record.split('\t')[3]

        if (gen1 == gen2):
            count += 1
            result_file.write(str(pos1.chrom) +  " " + str(pos1.pos) + " " +  gen1)
           
        if (not pos1.next()):
            break
        
    result_file.write("Result " + str(count))
    result_file.close()
    return count



main_file = sys.argv[1]
open_snp_dir = sys.argv[2]
result_dir = sys.argv[3]

result_list = []
for dir_name, dirnames, file_names in os.walk(sys.argv[2]):
    for file_name in file_names:
        try:
            print "Start processing ", file_name

            if (not file_name.endswith("23andme.txt")):
                print "Skipping"
                continue
            
            result_file = file(os.path.join(result_dir, file_name), "wr")
            dist = calc_inv_dist(main_file, os.path.join(dir_name, file_name), result_file)
            result_list +=[(dist, file_name)]
            print "Finish proccessing with result ", dist
        except Exception:
            print "Exception!!! Stop processing"
       
       

result_list.sort()

for i in range(0, 9):
    print result_list[i]

