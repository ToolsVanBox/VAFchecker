#!/usr/bin/python

#Import modules
import vcf as pyvcf
import pysam
import argparse
import multiprocessing as mp
import queue
import time
import sys
import collections
import subprocess

__version__ = subprocess.check_output(["git", "describe"]).strip().decode('UTF-8')

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-i', '--input', type=str, help='Input vcf file', required=True)
parser.add_argument('-b', '--bam', action='append', nargs="*", type=str, help='Input bam file', required=True)
parser.add_argument('-t', '--threads', default=4, type=int, help="Number of threads (default: %(default)s)")
parser.add_argument('-Q', '--QUAL', default=50, type=int, help="Report only variants with a minimal QUAL flag (default: %(default)s)")
parser.add_argument('-m', '--mapq', default=0, type=int, help="Include only reads with a minimal mapq (default: %(default)s)")
parser.add_argument('-p', '--base_phred_quality', default=0, type=int, help="Include only bases with a minimal base phred quality (default: %(default)s)")
parser.add_argument('-v', '--version', action='version', version=__version__)

args = parser.parse_args()
args.bam = [x for l in args.bam for x in l] #Flatten input list of bam files

def get_command_line():
    cmdline = sys.argv[0]
    for arg in vars(args):
        if type(getattr(args,arg)) is list:
            for a in getattr(args,arg):
                    cmdline += " --"+arg+" "+str(a)
        else:
            cmdline += " --"+arg+" "+str(getattr(args,arg))
    return( '"'+cmdline+'"' )

def fix_vcf_header( vcf_reader ):
    vcf_reader.infos['dbNSFP_clinvar_clnsig'] = pyvcf.parser._Info("dbNSFP_clinvar_clnsig",1,"String","Field 'clinvar_clnsig' from dbNSFP", None, None)
    return( vcf_reader )

def add_vcf_header( vcf_reader ):
    vcf_reader.formats['VAF'] = pyvcf.parser._Format('VAF',None,'Integer','Variant Allele Frequency calculated from the BAM file')
    vcf_reader.formats['CAD'] = pyvcf.parser._Format('CAD',None,'Integer','Calculated Allelic Depth, used for VAF calculation')
    vcf_reader.metadata['VAFcheckerCmd'] = [get_command_line()]
    return( vcf_reader )

def get_sample_name( bamfile ):
    header = bamfile.header
    if 'RG' in header:
        if type(header['RG']) is list:
            return(header['RG'][0]['SM'])
        else:
            return(header['RG']['SM'])
    return( False )

def update_call_data( call, edit_keys, edit_values ):
    f_keys = list(vcf_reader.formats.keys())
    d = dict(call.data._asdict())
    f_vals = list()
    for key in f_keys:
        if key in edit_keys:
            f_vals.append(edit_values[edit_keys.index(key)] )
        elif key in d:
            f_vals.append(d[key] )
        else:
            f_vals.append(None)
    handy_dict = dict(zip(f_keys, f_vals))
    f_keys.remove('GT')
    f_keys.insert(0,'GT')
    call.data = collections.namedtuple('CallData',f_keys)(**handy_dict)

def check_pileupread( pileupread ):
    if pileupread.alignment.is_duplicate:
        return( False )
    if pileupread.is_del:
        return( False )
    if pileupread.is_refskip:
        return( False )
    if not pileupread.query_position:
        return( False )
    if pileupread.alignment.mapq < args.mapq:
        return( False )
    if pileupread.alignment.query_qualities[pileupread.query_position] < args.base_phred_quality:
        return( False )

    return( True )

def check_record( record ):
    if record.QUAL < args.QUAL:
        return( False )
    if record.FILTER:
        return( False )

    return( True )

def parse_chr_vcf(q, q_out, chr_vcf_reader, bams):
    while True:
        try:
            contig = q.get(block=False,timeout=1)
            try:
                chr_vcf_reader.fetch(contig)
            except:
                continue
            for record in chr_vcf_reader.fetch(contig):
                if ( check_record( record ) ):
                    for call in (record.samples):
                        update_call_data(call, ['VAF','CAD'], [None, None])
                    for bam in bams:
                        F=pysam.AlignmentFile(bam,'rb')
                        sample_name = get_sample_name(F)
                        dv = [0]*len(record.ALT)
                        dr = 0
                        vaf = [0.0]*len(record.ALT)
                        for pileupcolumn in F.pileup(record.CHROM, int(record.POS)-1, int(record.POS), truncate=True, stepper='nofilter'):
                            for pileupread in pileupcolumn.pileups:
                                if ( check_pileupread( pileupread) ):
                                    for alt in record.ALT:
                                        if (len(record.REF) == 1 and len(alt) == 1):
                                            if pileupread.alignment.query_sequence[pileupread.query_position] == record.REF:
                                                dr+=1
                                            elif pileupread.alignment.query_sequence[pileupread.query_position] == alt:
                                                dv[record.ALT.index(alt)]+=1
                                        elif (len(record.REF) > 1 and len(alt) == 1):
                                            if ( pileupread.indel*-1 == len(record.REF)-1 ):
                                                dv[record.ALT.index(alt)]+=1
                                            elif pileupread.indel == 0:
                                                dr+=1
                                        elif ( len(record.REF) == 1 and len(alt) > 1 ):
                                            if ( pileupread.indel == len(alt)-1 ):
                                                dv[record.ALT.index(alt)]+=1
                                            elif pileupread.indel == 0:
                                                dr+=1
                                        else:
                                            if ( pileupread.indel == (len(alt)-len(record.REF)) ):
                                                dv[record.ALT.index(alt)]+=1
                                            elif pileupread.indel == 0:
                                                dr+=1
                        for x in range(0,len(dv)):
                            try:
                                vaf[x] = float("{0:.2f}".format(dv[x]/float(dv[x]+dr)))
                            except ZeroDivisionError:
                                continue
                        for call in (record.samples):
                            if call.sample == sample_name:
                                cad = list(dv)
                                cad.insert(0,dr)
                                update_call_data(call, ['VAF','CAD'], [vaf, cad])
                    format_list = list(vcf_reader.formats.keys())
                    format_list.remove('GT')
                    format_list.insert(0,'GT')
                    record.FORMAT = ":".join(format_list)
                    vcf_writer.write_record(record)
        except queue.Empty:
            break
    q_out.put('done')


vcf_reader = pyvcf.Reader(filename=args.input)
vcf_reader = fix_vcf_header(vcf_reader)

vcf_reader = add_vcf_header(vcf_reader)

vcf_writer = pyvcf.Writer(open('/dev/stdout', 'w'), vcf_reader)

contig_list = []
for contig in vcf_reader.contigs:
    contig_list.append(contig)

q = mp.Queue()
q_out = mp.Queue()
for contig in contig_list:
    q.put(contig)

processes = [mp.Process(target=parse_chr_vcf, args=(q, q_out, vcf_reader, args.bam)) for x in range(args.threads)]

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
