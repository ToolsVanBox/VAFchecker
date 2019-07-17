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
import os
import glob

# Get version from git
#__version__ = subprocess.check_output(["git", "describe"]).strip().decode('UTF-8')
__version__ = 'v1.1.2'

# Set arguments
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-i', '--input', type=str, help='Input vcf file', required=True)
parser.add_argument('-b', '--bam', action='append', nargs="*", type=str, help='Input bam file', required=True)
parser.add_argument('-t', '--threads', default=4, type=int, help="Number of threads (default: %(default)s)")
parser.add_argument('-c', '--clonal_threshold', default=0.3, type=float, help="Sample reported as subclonal if VAF is lower (default: %(default)s)")
parser.add_argument('-a', '--absent_threshold', default=0.0, type=float, help="Sample reported as absent if VAF is lower(default: %(default)s)")
parser.add_argument('-Q', '--QUAL', default=50, type=int, help="Report only variants with a minimal QUAL flag (default: %(default)s)")
parser.add_argument('-m', '--mapq', default=0, type=int, help="Include only reads with a minimal mapq (default: %(default)s)")
parser.add_argument('-p', '--base_phred_quality', default=0, type=int, help="Include only bases with a minimal base phred quality (default: %(default)s)")
parser.add_argument('-v', '--version', action='version', version=__version__)

args = parser.parse_args()
# Flatten input list of bam files
args.bam = [x for l in args.bam for x in l]

def get_command_line():
    """
    Function to get the commandline arguments
    Return: A string with the actual command.
    """
    cmdline = sys.argv[0]
    for arg in vars(args):
        if type(getattr(args,arg)) is list:
            for a in getattr(args,arg):
                    cmdline += " --"+arg+" "+str(a)
        else:
            cmdline += " --"+arg+" "+str(getattr(args,arg))
    return( '"'+cmdline+'"' )

def fix_vcf_header( vcf_reader ):
    """
    Function to fix fields in the vcf header
    Input: A vcf reader object
    Return: The vcf reader object with fixed headers
    """
    #dbNSFP_clinvar_clnsig has a Integer type but sometimes it is a String, e.g. 2|2
    vcf_reader.infos['dbNSFP_clinvar_clnsig'] = pyvcf.parser._Info("dbNSFP_clinvar_clnsig",1,"String","Field 'clinvar_clnsig' from dbNSFP", None, None)
    return( vcf_reader )

def add_vcf_header( vcf_reader ):
    """
    Function to add a new field to the vcf header
    Input: A vcf reader object
    Return: The vcf reader object with new headers added
    """
    vcf_reader.formats['VAF'] = pyvcf.parser._Format('VAF',None,'Integer','Variant Allele Frequency calculated from the BAM file')
    vcf_reader.formats['CAD'] = pyvcf.parser._Format('CAD',None,'Integer','Calculated Allelic Depth, used for VAF calculation')
    vcf_reader.metadata['VAFcheckerCmd'] = [get_command_line()]
    vcf_reader.infos['ABSENT'] = pyvcf.parser._Info('ABSENT',1,'Integer','Number of samples without the variant', None, None)
    vcf_reader.infos['SUBCLONAL'] = pyvcf.parser._Info('SUBCLONAL',1,'Integer','Number of samples with a subclonal variant', None, None)
    vcf_reader.infos['CLONAL'] = pyvcf.parser._Info('CLONAL',1,'Integer','Number of samples with a clonal variant', None, None)
    vcf_reader.infos['ABSENT_SAMPLES'] = pyvcf.parser._Info('ABSENT_SAMPLES',None,'String','Samples without the variant', None, None)
    vcf_reader.infos['SUBCLONAL_SAMPLES'] = pyvcf.parser._Info('SUBCLONAL_SAMPLES',None,'String','Samples with a subclonal variant', None, None)
    vcf_reader.infos['CLONAL_SAMPLES'] = pyvcf.parser._Info('CLONAL_SAMPLES',None,'String','Samples with a clonal variant', None, None)
    return( vcf_reader )

def get_sample_name( bamfile ):
    """
    Function to get the sample name from the bam file
    Input: An AlignmentFile object of the bam file
    Return: The sample or False if there is no SM tag in the bam header
    """
    header = bamfile.header
    if 'RG' in header:
        if type(header['RG']) is list:
            return(header['RG'][0]['SM'])
        else:
            return(header['RG']['SM'])
    return( False )

def update_call_data( call, edit_keys, edit_values ):
    """
    Function to add or update a field to the format field.
    This will be automatically update in the call object
    Input: A call object
    Input: A list with format fields
    Input: A list with format values
    """
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
    """
    Function to check a pileup read.
    Returns True if the read needs to be kept and returns False if read can be skipped.
    Input: Pileupread object
    Return: True or False
    """
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
    """
    Function to check a record.
    Returns True if the record needs to be kept and returns False if record can be skipped.
    Input: Record object
    Returns: True or False
    """
    if record.QUAL < args.QUAL:
        return( False )
    if record.FILTER:
        return( False )
    return( True )

def parse_chr_vcf(q, q_out, contig_vcf_reader, bams):
    """
    Function to parse the vcf per contig.
    Write the new record to a vcf file.
    Input: Queue object
    Input: Queue out object
    Input: VCF reader object
    Input: List with the bam names
    """
    while True:
        try:
            # Get contig one by one from the queue
            contig = q.get(block=False,timeout=1)
            contig_vcf_writer = pyvcf.Writer(open("./VAFchecker_tmp/"+contig+".vcf",'w'), contig_vcf_reader)
            try:
                # Try to parse the specific contig from the vcf
                contig_vcf_reader.fetch(contig)
            except:
                # Skip contig if this one is not present in the vcf file
                continue
            for record in contig_vcf_reader.fetch(contig):
                clonal_samples = [[]]*len(record.ALT)
                subclonal_samples = [[]]*len(record.ALT)
                absent_samples = [[]]*len(record.ALT)
                if ( check_record( record ) ):
                    for call in (record.samples):
                        # Add empty VAF and CAD tag to the record
                        update_call_data(call, ['VAF','CAD'], [None, None])
                    for bam in bams:
                        F=pysam.AlignmentFile(bam,'rb')
                        sample_name = get_sample_name(F)
                        dv = [0]*len(record.ALT)
                        dr = 0
                        vaf = [0.0]*len(record.ALT)
                        for pileupcolumn in F.pileup(record.CHROM, int(record.POS)-1, int(record.POS), truncate=True, stepper='nofilter',min_base_quality=args.base_phred_quality):
                            for pileupread in pileupcolumn.pileups:
                                if ( check_pileupread( pileupread) ):
                                    for alt in record.ALT:
                                        # If variant is a SNV
                                        if (len(record.REF) == 1 and len(alt) == 1):
                                            if pileupread.alignment.query_sequence[pileupread.query_position] == record.REF:
                                                dr+=1
                                            elif pileupread.alignment.query_sequence[pileupread.query_position] == alt:
                                                dv[record.ALT.index(alt)]+=1
                                        # If variant is an INDEL, in this case a deletion
                                        elif (len(record.REF) > 1 and len(alt) == 1):
                                            if ( pileupread.indel*-1 == len(record.REF)-1 ):
                                                dv[record.ALT.index(alt)]+=1
                                            elif pileupread.indel == 0:
                                                dr+=1
                                        # If variant is an INDEL, in this case an insertion
                                        elif ( len(record.REF) == 1 and len(alt) > 1 ):
                                            if ( pileupread.indel == len(alt)-1 ):
                                                dv[record.ALT.index(alt)]+=1
                                            elif pileupread.indel == 0:
                                                dr+=1
                                        # If variant is an INDEL
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
                                for vaf_idx in range(len(vaf)):
                                    if vaf[vaf_idx] <= args.absent_threshold:
                                        absent_samples[vaf_idx].append(call.sample)
                                    elif vaf[vaf_idx] < args.clonal_threshold:
                                        subclonal_samples[vaf_idx].append(call.sample)
                                    else:
                                        clonal_samples[vaf_idx].append(call.sample)
                    format_list = list(vcf_reader.formats.keys())
                    format_list.remove('GT')
                    format_list.insert(0,'GT')
                    record.FORMAT = ":".join(format_list)
                    record.INFO['ABSENT'] = [len(x) for x in absent_samples]
                    record.INFO['SUBCLONAL'] = [len(x) for x in subclonal_samples]
                    record.INFO['CLONAL'] = [len(x) for x in clonal_samples]
                    record.INFO['ABSENT_SAMPLES'] = ["|".join(x) for x in absent_samples]
                    record.INFO['SUBCLONAL_SAMPLES'] = ["|".join(x) for x in subclonal_samples]
                    record.INFO['CLONAL_SAMPLES'] = ["|".join(x) for x in clonal_samples]
                    contig_vcf_writer.write_record(record)
        # Break the loop if the queue is empty
        except queue.Empty:
            break
    q_out.put( 'done' )

# Read the vcf, fix and add fields to the header
vcf_reader = pyvcf.Reader(filename=args.input)
vcf_reader = fix_vcf_header(vcf_reader)
vcf_reader = add_vcf_header(vcf_reader)

# Open a vcf for writing
vcf_writer = pyvcf.Writer(open('/dev/stdout', 'w'), vcf_reader)

try:
    os.stat('./VAFchecker_tmp')
except:
    os.mkdir('./VAFchecker_tmp')


# Read the contig fields in vcf header to get all contigs
contig_list = []
for contig in vcf_reader.contigs:
    contig_list.append(contig)

# Create an input queue with the contigs and an empty output queue
q = mp.Queue()
q_out = mp.Queue()
for contig in contig_list:
    q.put(contig)

# Create number of processes to parse the vcf file
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
# Give tasks a chance to put more data in
    time.sleep(5)
    if not q.empty():
        continue
    liveprocs = [p for p in liveprocs if p.is_alive()]

for p in processes:
    p.join()

for tmp_vcf in glob.glob('./VAFchecker_tmp/*.vcf'):
    tmp_vcf_reader = pyvcf.Reader(filename=tmp_vcf)
    for record in tmp_vcf_reader:
        vcf_writer.write_record(record)
    os.remove(tmp_vcf)

#os.rmdir('./VAFchecker_tmp')
