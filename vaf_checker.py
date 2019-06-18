#!/usr/bin/python

#Import modules
import vcf as pyvcf
import pysam
import argparse
import multiprocessing as mp
import queue
import time

__version__ = '1.0.0'

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-i', '--input', type=str, help='Input vcf file', required=True)
parser.add_argument('-b', '--bam', action='append', nargs="*", type=str, help='Input bam file', required=True)
parser.add_argument('-t', '--threads', default=4, type=int, help='Number of threads [default: 4]')
parser.add_argument('-q', '--qual', default=50, type=int, help='Report only variants with a minimal qual [default: 50]')
#parser.add_argument('-f', '--filter', default='PASS', action='append', help='Report only variants with one of the filter [default: PASS]')
parser.add_argument('-m', '--mapq', default=0, type=int, help='Include only reads with a minimal mapq [default: 0]')
parser.add_argument('-p', '--base_phred_quality', default=37, type=int, help='Include only bases with a minimal base phred quality [default: 37]')
parser.add_argument('-v', '--version', action='version', version=__version__)

args = parser.parse_args()
args.bam = [x for l in args.bam for x in l] #Flatten input list of bam files

print("Variant"+"\t"+"\t".join(args.bam))

def parse_chr_vcf(q, q_out, vcf, bams):
    chr_vcf_reader = pyvcf.Reader(filename=vcf)
    chr_vcf_reader.infos['dbNSFP_clinvar_clnsig'] = pyvcf.parser._Info("dbNSFP_clinvar_clnsig",1,"String","Field 'clinvar_clnsig' from dbNSFP", None, None)

    while True:
        try:
            contig = q.get(block=False,timeout=1)
            try:
                chr_vcf_reader.fetch(contig)
            except:
                continue
            for record in chr_vcf_reader.fetch(contig):
                if len(record.ALT) > 1 or record.QUAL < args.qual or record.FILTER:
                    continue
                vafs = []
                for bam in bams:
                    F=pysam.AlignmentFile(bam,'rb')
                    dv = 0
                    dr = 0
                    for pileupcolumn in F.pileup(record.CHROM, int(record.POS)-1, int(record.POS), truncate=True, stepper='nofilter'):
                        for pileupread in pileupcolumn.pileups:
                            if pileupread.alignment.is_duplicate:
                                continue
                            if pileupread.query_position and pileupread.alignment.mapq >= args.mapq and pileupread.alignment.query_qualities[pileupread.query_position] >= args.base_phred_quality:
                                if pileupread.is_del:
                                    print("DEL")
                                if not pileupread.is_del and not pileupread.is_refskip:
                                    if (len(record.REF) == 1 and len(record.ALT[0]) == 1):
                                        if pileupread.alignment.query_sequence[pileupread.query_position] == record.REF:
                                            dr+=1
                                        elif pileupread.alignment.query_sequence[pileupread.query_position] in record.ALT:
                                            dv+=1
                                    else:
                                        if ( len(record.REF) > 1 ):
                                            if ( pileupread.indel*-1 == len(record.REF)-1 ):
                                                dv+=1
                                            elif pileupread.indel == 0:
                                                dr+=1
                                        elif ( len(record.ALT[0]) > 1 ):
                                            if ( pileupread.indel == len(record.ALT[0])-1 ):
                                                dv+=1
                                            elif pileupread.indel == 0:
                                                dr+=1
                    try:
                        vafs.append(str(dv/float(dv+dr)))
                    except ZeroDivisionError:
                        vafs.append('0.0')
                print( record.CHROM+"_"+str(record.POS)+"_"+record.REF+"/"+str(record.ALT[0])+"\t"+"\t".join(vafs) )




        except queue.Empty:
            break
    q_out.put('done')


vcf_reader = pyvcf.Reader(filename=args.input)
contig_list = []
for contig in vcf_reader.contigs:
    contig_list.append(contig)

q = mp.Queue()
q_out = mp.Queue()
for contig in contig_list:
    q.put(contig)

processes = [mp.Process(target=parse_chr_vcf, args=(q, q_out, args.input, args.bam)) for x in range(args.threads)]

for p in processes:
    p.start()

liveprocs = list(processes)
while liveprocs:
    time.sleep(5)
    try:
        while 1:
            done = q_out.get(block=False, timeout=1)
    except queue.Empty:
        pass

    time.sleep(5)    # Give tasks a chance to put more data in
    if not q.empty():
        continue
    liveprocs = [p for p in liveprocs if p.is_alive()]

for p in processes:
    p.join()
