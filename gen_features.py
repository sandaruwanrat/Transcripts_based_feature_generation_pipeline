import os
import sys
import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import numpy as np
from io import StringIO 
import pyranges as pr
from numpy import array
#import pybedtools
import time
from Bio.SeqUtils import gc_fraction
from itertools import product
from email.policy import default


####################################################################################################################

parser= argparse.ArgumentParser(description="Generate sequnce based features",
                                formatter_class=argparse.RawDescriptionHelpFormatter)



# parser.add_argument('-F','--fa_dir',default="")
# parser.add_argument('-S','--sample',default="")
# parser.add_argument('-T','--transgene',default="")

parser.add_argument('-F', '--fa_dir', required=True, help="Directory containing FASTA file")
parser.add_argument('-S', '--sample', required=True, help="Sample FASTA file name")
parser.add_argument('-T', '--transgene', required=True, help="Transgene name or ID")

args=parser.parse_args()


round_dir= args.fa_dir
fa_sample= args.sample
tg= args.transgene



#count di neucleotide count
def Counting_di(seq):
    sequence = seq
    di_nucleotides = [''.join(x) for x in product('ACGT', repeat=2)]
    # Generate all possible di-nucleotides
    di_count = {di: 0 for di in di_nucleotides} # Initialize dictionary with all di-nucleotides
  
    for i in range(len(sequence) - 1):
        di = sequence[i:i+2].upper() # Extract current di-nucleotide
        if di in di_count:
            di_count[di] += 1 # Increment count if di-nucleotide is valid
    
    #print(di_count)

    return di_count



#count tri neucleotide count
def Counting_tri(seq):
    sequence = seq
    tri_nucleotides = [''.join(x) for x in product('ACGT', repeat=3)] # Generate all possible tri-nucleotides
    tri_count = {tri: 0 for tri in tri_nucleotides} # Initialize dictionary with all tri-nucleotides
    
    for i in range(len(sequence) - 2):
        tri = sequence[i:i+3].upper() # Extract current tri-nucleotide
        if tri in tri_count:
            tri_count[tri] += 1 # Increment count if tri-nucleotide is valid
        

    return tri_count


#count tetra neucleotide count
def Counting_tetra(seq):
    sequence = seq
    tetra_nucleotides = [''.join(x) for x in product('ACGT', repeat=4)] # Generate all possible terai-nucleotides
    tetra_count = {tetra: 0 for tetra in tetra_nucleotides} # Initialize dictionary with all tri-nucleotides
    
    for i in range(len(sequence) - 3):
        tetra = sequence[i:i+4].upper() # Extract current tri-nucleotide
        if tetra in tetra_count :
            tetra_count[tetra] += 1 # Increment count if tri-nucleotide is valid
        

    return tetra_count




####################################################################################################################

       
in_fasta=str(fa_sample)+ '.'+ str(tg)+'.fasta'
print(in_fasta)
infa_pth=os.path.join(round_dir ,in_fasta)
print(infa_pth)
outdr = os.path.split(infa_pth)[0]
output_file = str(fa_sample)+ '_'+ str(tg) + '_3p5ptr_features_out.txt'
output_pth = os.path.join(outdr, output_file)
print(output_pth)




def run_fgen(input_dir,sname,infa_pth):
    first_iteration = True
    for seq_record in SeqIO.parse(infa_pth, 'fasta'):

        readid=seq_record.id
        #print(readid)
        sequence = str(seq_record.seq)
        gc_content=gc_fraction(sequence)
        
        read_dic= {'seq_id':seq_record.id,'length': len(sequence) , 'gc_cont': gc_content }
        di_dic=Counting_di(sequence)
        
        read_dic.update(di_dic)

        tri_dic=Counting_tri(sequence)
        
        read_dic.update(tri_dic)

        tetra_dic=Counting_tetra(sequence)
        
        read_dic.update(tetra_dic)


        ddf=pd.DataFrame.from_dict(read_dic,orient='index',columns=['Count'])

        dfn = (ddf.T).reset_index(drop=True)




        dfn.to_csv(output_pth, mode='a', header=first_iteration, index=False,sep='\t')
        if first_iteration:
            first_iteration=False
        




def main():
    start= time.time()
    run_fgen(input_dir=round_dir,sname=fa_sample,infa_pth=infa_pth)
    etime=(time.time()-start)
    
    
    
    

    

if __name__ == "__main__":
    main()





