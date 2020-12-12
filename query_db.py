from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse
import numpy as np
import dna_jellyfish as jellyfish



def sequence_parser(database, jellypath, kmerSize, output):
    jelly = jellyfish.QueryMerFile(jellypath)
    out = open(output, 'w')
    with open(database, 'r') as db:
        for l in db:
            sequence = l.split()[0]

            for i in range(len(sequence) - kmerSize + 1):  #split the sequence of the database into 51bp kmers size and move one position each time.  
                kmerF = sequence[i:i+kmerSize].upper() #put all the K-mers from i:i with the size and upper in the variable KmerF  
                mer = jellyfish.MerDNA(kmerF) # Now read the sequence in the KmerF variable that we just created using this MerDNA#https://github.com/gmarcais/Jellyfish/blob/master/swig/Readme.md
                mer.canonicalize() #keep the representative kmers in kmerF
                if jelly[mer]: #if the K-mer in the (jelly variable) is in the sequence haplotype DB (mer), will write the complete line (l) because this shlud match complete in the output?.
                    out.write(l) #It writes the whole line of the haplotype since here there is no indications to use and analyze the other columns in the DB file. Therefore if there is just one K-mer present, does this mean that all the haplotype is present? or just partially? if so, are we capturing false positive K-mers-haplotypes in the Vectors,
                    break
    
    out.close()
       


if __name__ == '__main__': # when python interpreter read a file sets a variable 'name' and then executes the code. To be read in the comand line lik python query_db.py. This happen becaus this modfule is the main program

    parser = argparse.ArgumentParser(description = "Parse haplotype database.")

    db_params = parser.add_argument_group('Parameters of database to be parsed.')
    db_params.add_argument('-d', '--database',  required=True, help='Path to database.')

    jelly_params = parser.add_argument_group('Parameters of query jellyfish.')
    jelly_params.add_argument('-k', '--kmersize',  required=True, help='Kmer size specified while making jellyfish.')
    jelly_params.add_argument('-j', '--jellypath',  required=True, help='Path to jellyfish of query file.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', help='Path of output file.')

    args = parser.parse_args()

    sequence_parser(args.database, args.jellypath, int(args.kmersize), args.output)