#!/bin/bash
"exec" "`dirname $0`/../python_env/bin/python" "$0" "$@"

# ^^^
# the cmd above ensures that the correct python environment is 
# selected to execute this script.
# The correct environment is the one belonging to uap, since all 
# neccessary python modules are installed there.



#post_cufflinks_merge.py

# test it with
# $ tools/post_cufflinks_merge.py /data/bioinf/projects/data/2015_mouseTcells_BBP/cuffcompare_segemehl_local/magic.combined.gtf 


import sys
import re
import argparse
import pprint
#import yaml
from collections import OrderedDict 
from collections import defaultdict
import numpy
from Bio import SeqIO

pp = pprint.PrettyPrinter(indent=4)
  


def read_arguments():
    parser = argparse.ArgumentParser(description="Filters merged.gtf from cufflinks for novel transcripts and spits out some statistics")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin, help ="Infile: default=stdin")

    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                    default=sys.stdout, help ="Outfile: default=stdout")

    parser.add_argument('--logfile', nargs='?', type=argparse.FileType('w'),
                    default=sys.stderr, help ="Novel transcript  metrics file:  default=stderr ")
        
    
    parser.add_argument('--remove-unstranded',  action='store_true',
                    default=False, help ="Removes transcripts without strand specifity")

    parser.add_argument('--remove-gencode',  action='store_true',
                    default=False, help ="hard removal of gtf line which match 'ENS' in gene_name field")

    parser.add_argument('--string',  
                        default='ENS', help ="string to match in gtf field gene_name for discarding")


    parser.add_argument('--remove-by-gene-name',  action='store_true',
                        default=False, help ="remove gtf if matches string in gene_name field")

    parser.add_argument('--filter-by-class',  action='store_true',
                        default=False, help ="remove gtf if any class is found in class_code field, requieres class_list")
 
    parser.add_argument('--class-list',  
                        default=None, help ="class codes to be removed possible '=,c,j,e,i,o,p,r,u,x,s,.'")

    parser.add_argument('--filter-by-class-and-gene-name',  action='store_true',
                        default=False, help ="combines remove-by-class and remove-by-gene-name")

    return parser.parse_args()



def eval_arguments(args):
    if args.class_list:
        t = args.class_list
        t = t.replace("'", '')
        args.class_list = t.split(',')
    return(args)




## called from main                     

def output_metrics (metrics):
#    args.logfile.write( yaml.dump(metrics, default_flow_style=False))
    pass




def write_t_obj(t_obj):
    gtf_dicts = t_obj['gtf'] 

    for i in gtf_dicts:
        t_list=[]
        keys = i.keys()
        for key in keys[8:]:
            opt_id = key
            opt_val = str(i[key])
            merge = opt_id + ' ' + '"' + opt_val + '"'
            t_list.append(merge)

        line_end = ';'.join(t_list)
        
        t_list=[]
        for key in keys[:8]:
            opt_val = str(i[key])
            t_list.append(opt_val)

        t_list.append(line_end)
        line = '\t'.join(t_list)
        args.outfile.write(line + '\n')


def read_gtf_line(line):
    """
    Returns the gtf/gff line in a dictionary.
    hickups if no 'gene_id', 'transcript_id' are present.
    
    http://genome.ucsc.edu/FAQ/FAQformat.html#format3

    seqname - The name of the sequence. Must be a chromosome or scaffold.
    source - The program that generated this feature.
    feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    start - The starting position of the feature in the sequence. The first base is numbered 1.
    end - The ending position of the feature (inclusive).
    score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
    strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
    group - A...

    gene_id value - A globally unique identifier for the genomic source of the sequence.
    transcript_id value

    
    chr1 Cufflinks exon 798899 799191 . + . gene_id "XLOC_000024"; transcript_id "TCONS_00000072"; exon_number "1"; oId "CUFF.8.1"; class_code "u"; tss_id "TSS51";

    """

    fields = line.split('\t')
    keylist = 'seqname, source, feature, start, end, score, strand, frame'.split(', ')

    optional_field = fields[8]
    gtf_dict = OrderedDict(zip(keylist, fields[0:8]))


        
    opt_entries = optional_field.split(';')


    for opt_entry in opt_entries:
        if opt_entry == '':
            continue
        else:
            opt_id, opt_val = str.split(opt_entry)[0:2]
            gtf_dict[opt_id] = opt_val.replace('"', '')

    
    for key in ['gene_id', 'transcript_id']:
        if not key in gtf_dict:
            raise Exception('key:{0} not in gtf_line\n {1}'.format(key, line))
            
    for numbers in ['start', 'end', 'exon_number']:
        gtf_dict[numbers] = int(gtf_dict[numbers])

        
    return gtf_dict
    

def make_transcript_object(gtf_dicts):

    t_obj =  defaultdict(lambda: defaultdict(dict))
    t_obj['info']['n_exon'] = 0
    t_obj['info']['class_code'] = ()
    t_obj['info']['exon_lengths'] = []
    t_obj['info']['transcript_length'] = 0
    t_obj['info']['seqname'] = ''
    t_obj['info']['strand'] = ''
    t_obj['info']['ok'] = True
    t_obj['info']['discard'] = False

    t_obj['gtf'] = []

    
    chrom       = []
    strands     = []
    classes     = []

    for gtf_dict in gtf_dicts:
        t_obj['info']['n_exon'] += 1
        t_obj['info']['exon_lengths'].append( gtf_dict['end'] - gtf_dict['start'] )

        chrom.append(gtf_dict['seqname'])
        strands.append(gtf_dict['strand'])
        classes.append(gtf_dict['class_code'])
        
        t_obj['gtf'].append(gtf_dict)

        

    t_obj['info']['transcript_length'] = sum(t_obj['info']['exon_lengths'])
    
    
    if not all_same(strands):
        info['discard'] = True
    else:
        t_obj['info']['strand'] = strands[0]
        
    if not all_same(chrom):
        info['discard'] = True
    else:
        t_obj['info']['seqname'] = chrom[0]

    if not all_same(classes):
        t_obj['info']['class_code'] = 'mixed'
    else:
        t_obj['info']['class_code'] = classes[0]
        
    
    return t_obj

def collect(gtf_dicts, gtf_dict, ID):
    if ID == None:
        ID = gtf_dict['transcript_id']
        
    if ID ==  gtf_dict['transcript_id']:
        gtf_dicts.append(gtf_dict)
        process_switch = 0
    else: 
        process_switch = 1
        ID = gtf_dict['transcript_id']
     
    return gtf_dicts, gtf_dict, ID, process_switch


def remove_unstranded(gtf_dict):
    discard = False
    
    if gtf_dict['strand'] == '.':
        discard = True

    return discard 
      


def add_metrics(metrics, t_obj):
    metrics['transcripts'] += 1

    if t_obj['info']['discard']:
        metrics['discard'] += 1
        return metrics

    tmp_val = t_obj['info']['seqname'] 
    metrics['transcripts_chr'][tmp_val] += 1 

    tmp_val = t_obj['info']['strand'] 
    metrics['strands'][tmp_val] += 1
   
    tmp_val = t_obj['info']['class_code'] 
    metrics['class_codes'][tmp_val] += 1

    tmp_val = t_obj['info']['n_exon'] 
    metrics['exons'][tmp_val] += 1


    tmp_val = t_obj['info']['transcript_length'] 
    metrics['transcript_length'].append(tmp_val)

    tmp_list = t_obj['info']['exon_lengths'] 
    metrics['exon_length'].extend(tmp_list)
    

    return metrics

def init_metrics():
    pass
    """Returns a dict like thingy for counting"""
    #new test add 

    metrics =  defaultdict(lambda: defaultdict(dict))

    metrics['transcripts'] = 0
    metrics['discard'] = 0
    metrics['transcripts_chr'] = defaultdict(int)
    metrics['exons'] = defaultdict(int)
    metrics['strands'] =  defaultdict(int)
    metrics['class_codes'] =  defaultdict(int)
    metrics['exon_length'] = []
    metrics['transcript_length'] = []



    return metrics
     
def get_averages_metrics(metrics):

    mean_exon_length         = int(numpy.mean(metrics['exon_length']))
    median_exon_length       = int(numpy.median(metrics['exon_length']))
    mean_transcript_length   = int(numpy.mean(metrics['transcript_length']))
    median_transcript_length = int(numpy.median(metrics['transcript_length']))
    

    metrics['mean_exon_length']            = mean_exon_length
    metrics['median_exon_length']          = median_exon_length
    metrics['mean_transcript_length']      = mean_transcript_length
    metrics['median_transcript_length']      = median_transcript_length

    metrics.pop("exon_length", None)
    metrics.pop("transcript_length", None)


    return metrics


def all_same(items):
    """
    tests if all items in alist are the same
    http://stackoverflow.com/a/3787983
    """ 
    return all(x == items[0] for x in items)



def remove_gencode(gtf_dict):
    pattern = re.compile('ENS')
    discard = False

    if 'gene_name' in gtf_dict:
        result = pattern.search(gtf_dict['gene_name'])
        if result:
            discard = True

    return discard 




def remove_by_gene_name(gtf_dict, string):
    pattern = re.compile(string)
    discard = False

    if 'gene_name' in gtf_dict:
        result = pattern.search(gtf_dict['gene_name'])
        if result:
            discard = True

    return discard 



def filter_by_class(gtf_dict, remove_class_list):
    """ from tophat
    Priority Code Description
    1 = Complete match of intron chain
    2 c Contained
    3 j Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript
    4 e Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment.
    5 i A transfrag falling entirely within a reference intron
    6 o Generic exonic overlap with a reference transcript
    7 p Possible polymerase run-on fragment (within 2Kbases of a reference transcript)
    8 r Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case
    9 u Unknown, intergenic transcript
    10 x Exonic overlap with reference on the opposite strand
    11 s An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)
    12 . (.tracking file only, indicates multiple classifications)"""

    discard = False
    if 'class_code' in gtf_dict:
        if gtf_dict['class_code'] in remove_class_list:
            discard = True
            return discard 
    
def filter_by_class_and_gene_name(gtf_dict, string=None, class_list=None):
    discard = False
    gene_res  = remove_by_gene_name(gtf_dict, string)
    class_res = filter_by_class(gtf_dict, class_list)

    if (gene_res == True and class_res == True):
        discard = True
        return discard 
    

def main(args):
    ID = None
    # print(args)
    out = args.outfile
    gtf_dicts =[]
    metrics = init_metrics()

    for line in args.infile:
        if line[0] == '#':
            continue
        line = line.rstrip('\n')
        gtf_dict = read_gtf_line(line)

                               

        if args.remove_gencode:
            if remove_gencode(gtf_dict):
                continue

        if args.remove_by_gene_name:
            if remove_by_gene_name(gtf_dict, args.string):
                continue 

        if args.filter_by_class:
            if filter_by_class(gtf_dict, args.class_list):
                continue 

        if args.filter_by_class_and_gene_name:
            if filter_by_class_and_gene_name(gtf_dict, args.string, args.class_list):
                continue

        if args.remove_unstranded:
            if remove_unstranded(gtf_dict):
                continue 

            

        #take all gtf entries belonging to  same 'transcript_id'
        (gtf_dicts, gtf_dict, ID, process_switch) =  collect(gtf_dicts, gtf_dict, ID)
        
                
        #new transcript encountered start processesing  
        if process_switch ==  1 :
            t_obj = make_transcript_object(gtf_dicts)
            metrics = add_metrics(metrics, t_obj)
            write_t_obj(t_obj)
            gtf_dicts =[]
            gtf_dicts.append(gtf_dict)

         

    #last block after while iteration        
    if not gtf_dicts:
        raise Exception("You discarded everthing") 
    t_obj = make_transcript_object(gtf_dicts                                   )
    metrics = add_metrics(metrics, t_obj)
    write_t_obj(t_obj)
    metrics = get_averages_metrics(metrics)
    output_metrics(metrics)

if __name__ == '__main__':
    args = read_arguments()
    args = eval_arguments(args)
    main(args)

__version__ = '0.001'
