from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse
import numpy as np
import dna_jellyfish as jellyfish
# import resource

def get_kmer_presence(kmerF):
    mer = jellyfish.MerDNA(kmerF)
    mer.canonicalize() 
    kmer_pres = []
    pres = int(qjellies[mer]>0)
    if pres:
        return None
    else:
        kmer_pres.append(1)
    return kmer_pres
    

def sequence_parser(assembly, kmerSize, output):
    out = open(output, 'w')
    fasta_sequences = SeqIO.parse(open(assembly),'fasta')
    for fasta in fasta_sequences:
        seq_name, sequence = fasta.id, str(fasta.seq)
        seq_split = seq_name.split(':')
        name = seq_split[0]
        sort_index = int(seq_split[1])
        
        haplotype = None
        for i in range(len(sequence) - kmerSize + 1):
            kmerF = sequence[i:i+kmerSize].upper()
            if 'N' in kmerF:
                continue
            kmer_pres = get_kmer_presence(kmerF)
            # print(kmer_pres)

            if haplotype is not None: 
                if kmer_pres == kmer_pres_prev:
                    haplotype += kmerF[-1]
                else:
                    kmer_pres_prev_str=[str(pres) for pres in kmer_pres_prev]
                    out.write(haplotype+'\t'+str(len(haplotype))+'\t'+name+'\t'+str(sort_index)+'\t'+'\t'.join(kmer_pres_prev_str)+'\n')
                    # out.write(str(len(haplotype))+'\t'+str(sort_index)+'\t'+'\n') # with this reduce the size of the output file (from 1000 to 68bytes). only size and position os printed, but name must be useful
                    haplotype = None
                    #reset position of the fasta
                    sort_index = int(seq_split[1])

            if kmer_pres is not None and haplotype is None:
                haplotype = kmerF
                sort_index += i # move one position in the fasta seq
                kmer_pres_prev = kmer_pres # changed to [1] pattern
    out.close()
       
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Parse fasta file to check presence/absence of k-mers in given lines.")

    sequence_params = parser.add_argument_group('Parameters of file to be parsed.')
    sequence_params.add_argument('-a', '--assembly',  required=True, help='Fasta file to parse.')
    sequence_params.add_argument('-l', '--linename',  required=True, help='Name of the line containing the sequence.')

    jelly_params = parser.add_argument_group('Parameters of jellyfish files used for checking presence/absence.')
    jelly_params.add_argument('-k', '--kmersize',  required=True, help='Kmer size specified while making jellyfish.')
    jelly_params.add_argument('-j', '--qjellies',  required=True, help='kmers of the query sample.')
    
    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', help='Path of output file.')

    args = parser.parse_args()
    qjellies = jellyfish.QueryMerFile(args.qjellies)

    sequence_parser(args.assembly, int(args.kmersize), args.output)