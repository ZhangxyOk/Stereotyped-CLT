import regex
import sys
import argparse
import os
from collections import defaultdict
## import self define modules
sys.path.append("/mnt/data/home/lzz/project/Reconstruct_lineage/scripts/NewReTree")
from modules.SequenceParser import convert2string, fqTOfa, formFASTA, seqco, seqre, formFASTQ


# sample barcodes
BC_sample_dict = {
    "CACTCGACTCTCGCGT": "CBRAd5",
    "TCTCGTCGCAGTCTCT": "CBRAd5",
    "TCTGTATCTCTATGTG": "CBRAd6",
    "GCGCGCGCACTCTCTG": "CBRAd6",
    "ACAGTCGAGCGCTGCG": "CBRAd4",
    "GCTGAGACGACGCGCG": "CBRAd4"
}

def demulti_sample(sampleBC_dict, sequence, mismatch=2):
    found_sample = ""
    for BC, sampleName in sampleBC_dict.items():
        f_search_pattern = "(?e)(?P<Primer>%s){e<=%d}" % (BC, mismatch)
        r_search_pattern = "(?e)(?P<Primer>%s){e<=%d}" % (seqco(seqre(BC)), mismatch)
        if regex.search(f_search_pattern, sequence) or regex.search(r_search_pattern, sequence):
            found_sample = sampleName
            break
    return found_sample


def Parsers():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=str, nargs='?', required=True, help="Input the FASTQ format file")
    parser.add_argument('--label', type=str, nargs='?', required=True, help="Output file label. e.g: R1.clean, R2.clean, other.clean")
    parser.add_argument('--outDir', type=str, nargs='?', required=True, help="Output results directory")
    parser.add_argument('--mismatch', type=int, nargs='?', required=True, default=2, help="Primers mismatch permit")
    #parser.add_argument('--numcpu', type=int, nargs='?', required=True, help="Number of cpu needed")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    options = Parsers()
    # creat file handle
    sampleName_filehandle = {}
    for name in set(BC_sample_dict.values()):
        sampleName_filehandle[name] = open(os.path.join(options.outDir, name + "_" + options.label + ".fastq"), 'w')
    # unrecognizable file
    unrec_file = open(os.path.join(options.outDir, "unrecognizable_reads" + "_" + options.label + ".fastq"), "w")
    # start to demultiplexing
    inputFile = open(options.infile, 'r')
    for record in formFASTQ(inputFile):
        _, seq, _, _ = record
        de_sampleName = demulti_sample(BC_sample_dict, seq, options.mismatch)
        #print(de_sampleName)
        if de_sampleName:
            sampleName_filehandle[de_sampleName].write("{}\n{}\n{}\n{}\n".format(*record))
            sampleName_filehandle[de_sampleName].flush()
        else:
            unrec_file.write("{}\n{}\n{}\n{}\n".format(*record))
            unrec_file.flush()
    # close files
    for file_handle in sampleName_filehandle.values():
        file_handle.close()
    unrec_file.close()
    inputFile.close()