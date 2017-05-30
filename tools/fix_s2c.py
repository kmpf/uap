#!/bin/bash
"exec" "`dirname $0`/../python_env/bin/python" "$0" "$@"

# ^^^
# the cmd above ensures that the correct python environment is 
# selected to execute this script.
# The correct environment is the one belonging to uap, since all 
# neccessary python modules are installed there.



import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import sys
import argparse
import pprint
from collections import deque
from collections import Counter

pp = pprint.PrettyPrinter(indent=4)



def read_arguments():
    parser = argparse.ArgumentParser(description="Reads s2c sam (converter from segemehl) and repairs some entries to pass picard validate sam")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                    default=sys.stdin, help ="Infile default reads from stdin")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                    default=sys.stdout, help ="Outfile default writes to stdout")
    return parser.parse_args()


def read_samline(line):
    fields = line.split('\t')
    keylist = 'qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual'.split(', ')

    optional_fields = fields[ 11: ]
    samdict = dict(zip(keylist,fields[0:11]))
    samdict['flag'] = int(samdict['flag'])


    for optEntry in optional_fields:
        (optID, optDEF, optVal) = optEntry.split(':')
        samdict[optID] = [optVal, optDEF]
     
    return samdict


def make_samline(samdict):
    #should have impleneted opt fields in sep sub hash  opt future enhancement
    #this would spare me the fucking iterations
    keylist = 'qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual'.split(', ')

    opt_keys = []
    for key in samdict:
        if key not in keylist:
            opt_keys.append(key)

    opt_keys.sort()        

    pre_line =[]
    for i in keylist:
        pre_line.append(samdict[i])
        
    for i in opt_keys:
        opt = ':'.join([i,samdict[i][1],samdict[i][0]])
        pre_line.append(opt)

    
    samline= "\t".join(str(v) for v in pre_line)

    
    return samline
    


def collect(sam_hits, samdict,ID):
    if ID == None:
        ID = samdict['qname']
        
    if ID ==  samdict['qname']:
        sam_hits.append(samdict)
        process_switch = 0
    else: 
        process_switch = 1
        ID = samdict['qname']
     
    return sam_hits, samdict, ID, process_switch
     



def pre_process_sam_hits(sam_hits):
    #see if paired end or single and no mixture, multiple 
    
    
    for samdict in sam_hits:
        if 'XI' not in samdict:
            sys.stderr.write("Samline entry without XI entry")
            pp.pprint(samdict)
            exit(1)
            
         #sort all segemehl entries by XI
        sam_hits = sorted(sam_hits, key=lambda k: int(k['XI'][0])) 

        temporary_list =[]
        for samdict in sam_hits:
            temporary_list.append(samdict['XI'][0])
        counter_list = Counter(temporary_list)

          
        info   = dict()
        info['counter_list'] = counter_list
        info['single'], info['paired'], info['NH']  = 0,  0, 0
        

        for XI, how_often in sorted(counter_list.items()):
            info['NH'] += 1
            if (how_often == 1):
                info['single'] +=1
            elif (how_often == 2):
                info['paired'] +=1

        if info['single'] == 0 and info['paired'] == 0:
             sys.stderr.write("empty list something went wrong dying")
             exit(1)
        elif info['single'] > 0 and info['paired'] > 0:
            pp.pprint(sam_hits)
            pp.pprint(info)
            sys.stderr.write("mixed reads of type info['single'] and info['paired'] wtf dying srsly")
            exit(1)
        elif info['single'] > 0 and info['paired'] == 0:
            info['type'] = 'single'
        elif info['single'] == 0 and info['paired'] > 0:
            info['type'] = 'paired'


        return info, sam_hits





def process_sam_hits(info, sam_hits):
    return_list = []
    
    if info['type'] == 'single':
        for samdict in sam_hits:
            #set clear and set multiple flag ### for now just removing 
            if info['NH'] > 1:
               "temp"
               # samdict['flag'] = setBit(samdict['flag'],8)
            else: 
                samdict['flag'] = clearBit(samdict['flag'],8)
            #if  paired set mate unmapped 
            samdict['flag'] = test_set_flag(samdict['flag'], 0, 0, 3)
            return_list.append(samdict)


    elif info['type'] == 'paired':
        sam_hits = deque(sam_hits)
        #pp.pprint ('#####')
        #pp.pprint (info)

        for i in range (info['NH']):
            samdict_A = sam_hits.popleft()
            samdict_B = sam_hits.popleft()

            if samdict_A['XI'][0] !=  samdict_B['XI'][0]:
                sys.stderr.write("XI fields mixed up forgot sorting ??")
                exit(1)
            if info['NH'] > 1:
               "temp"
               # samdict_A['flag'] = setBit(samdict_A['flag'],8)
               # samdict_B['flag'] = setBit(samdict_B['flag'],8)
            else: 
                samdict_A['flag'] = clearBit(samdict_A['flag'],8)
                samdict_B['flag'] = clearBit(samdict_B['flag'],8)

            #set mate info    

                #not fail safe assuming there are no unmapped ones
                #fix later should not happen in segemehl
            if (samdict_A['rname'] == samdict_B['rname']):
                samdict_A['mrnm'] = '='
                samdict_B['mrnm'] = '='
            else:
                samdict_A['mrnm'] = samdict_B['rname']
                samdict_B['mrnm'] = samdict_A['rname']

            samdict_A['mpos'] = samdict_B['pos']
            samdict_B['mpos'] = samdict_A['pos']
            

            #if mate reverse set mat reverse flag in my info 
            if  testBit(samdict_B['flag'],4) > 0:
                samdict_A['flag'] = setBit(samdict_A['flag'],5)

            if  testBit(samdict_A['flag'],4) > 0:
                samdict_B['flag'] = setBit(samdict_B['flag'],5)


            return_list.append(samdict_A)
            return_list.append(samdict_B)
    else:
        sys.stderr.write("check check arrrgh")
        exit(1)

    return return_list






def test_set_flag(flag, bit_to_check, response, bit_to_set):
    if testBit(flag, bit_to_check) > response:
        flag = setBit(flag,bit_to_set) 
    return flag
                        
  

# testBit() returns a nonzero result, 2**offset, if the bit at 'offset' is one.

def testBit(int_type, offset):
    mask = 1 << offset
    return(int_type & mask)

# setBit() returns an integer with the bit at 'offset' set to 1.

def setBit(int_type, offset):
    mask = 1 << offset
    return(int_type | mask)

# clearBit() returns an integer with the bit at 'offset' cleared.

def clearBit(int_type, offset):
    mask = ~(1 << offset)
    return(int_type & mask)

# toggleBit() returns an integer with the bit at 'offset' inverted, 0 -> 1 and 1 -> 0.

def toggleBit(int_type, offset):
    mask = 1 << offset
    return(int_type ^ mask)

 
            
def output(sam_hits, outfile):
    for samdict in sam_hits:
        samline = make_samline(samdict)
        outline = ''.join([samline, "\n"])
        outfile.write(outline)



def main(args):
    ID = None
    # print(args)
    out = args.outfile
    sam_hits =[]
    
    for line in args.infile:
        if line[0] == '@':
            args.outfile.write(line)
            continue
        line = line.rstrip('\n')
        #args.outfile.write(line)
        samdict = read_samline(line)

        # take all reads same name
        (sam_hits, samdict, ID, process_switch) =  collect(sam_hits, samdict, ID)
        
        #new readname  encountered start processeing  
        if process_switch ==  1 :
            #look at XI paires and determine if only paired or single are in the collection
            info, sam_hits = pre_process_sam_hits(sam_hits)
            #
            sam_hits = process_sam_hits(info, sam_hits)
            output(sam_hits, args.outfile);
            #
            sam_hits = []
            sam_hits.append(samdict)


    #last block after while iteration        
    info, sam_hits = pre_process_sam_hits(sam_hits)
    sam_hits = process_sam_hits(info, sam_hits)
    output(sam_hits, args.outfile);



if __name__ == '__main__':
    args = read_arguments()
    main(args)
