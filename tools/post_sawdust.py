#!/usr/bin/env python
# post_sawdust.py

import sys
import re
import argparse
from argparse import RawTextHelpFormatter
import pprint
import yaml
from collections import OrderedDict
import numpy
from Bio import SeqIO

pp = pprint.PrettyPrinter(indent=4)


def read_arguments():
    parser = argparse.ArgumentParser(
        description="takes segemehl SAM output and merges split reads into one SAM line in addition some entries are tinkered with")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin, help="Infile: default=stdin")

    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help="Outfile: default=stdout")

    parser.add_argument(
        '--logfile',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stderr,
        help="Read count metrics file:  default=stderr ")

    parser.add_argument(
        '--outfile-splits',
        nargs='?',
        type=argparse.FileType('w'),
        default=None,
        help="Write not joinable split reads(i.e wrong chr etc)  in a separate file")

    parser.add_argument(
        '--distance',
        type=int,
        default=200000,
        help="Maximum distance between two fragments in a split read")

    parser.add_argument(
        '--filter-snp-calling',
        action='store_const',
        const='tunred_on',
        default=None,
        help="Filters out reads for snp calling, no multiple, no diff chr")

    parser.add_argument(
        '--genome',
        default=None,
        help="Fasta file, like: hg19.fa Only necessary for unstranded RNA Libraries")

    parser.add_argument(
        '--read-type',
        choices=[
            'single',
            'paired'],
        default='paired',
        help="single or paired end default=paired")

    parser.add_argument(
        '--library-type',
        choices=[
            'fr-unstranded',
            'fr-firststrand',
            'fr-secondstrand'],
        default=None,
        help="See tophat manual.")

    parser.add_argument(
        '--seq-type',
        choices=[
            'RNA',
            'DNA'],
        default='RNA',
        help="Specifies the molecule you sequenced commonly referred to as RNA- or DNA-seq")

    return parser.parse_args()


def eval_arguments(args):

    if args.seq_type is None:
        pass
        raise Exception(
            'Argument --seq-type is {0}; argument is required; choices: DNA, RNA'.format(args.seq_type))

    if args.seq_type == 'RNA':
        pass
        if args.library_type is None:
            my_choices = 'fr-unstranded', 'fr-firststrand', 'fr-secondstrand'
            raise Exception(
                'Argument --library-type is {0}; argument is required; choices: {1}'.format(
                    args.library_type, my_choices))

    if (args.library_type == 'fr-unstranded' and args.seq_type == 'RNA'):
        if args.genome is None:
            raise Exception(
                'Argument --genome is {0}; please specify fasta file '.format(args.genome))
        handle = open(args.genome, "rU")
        args.genome_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()


# called from main
def init_metrics():
    """Returns a dict like thingy for counting"""
    # new test add
    metrics = []

    indict = OrderedDict()
    split_types = [
        'not_same_chr',
        'not_same_strand',
        'max_dist_breached',
        'fragment_overlap',
        'order_switch',
        'OK']

    indict['templates'] = 0
    indict['unmapped'] = 0
    indict['multiple'] = 0
    indict['sub_multiple'] = {}
    indict['sub_multiple']['r1'] = 0
    indict['sub_multiple']['r2'] = 0

    indict['unique'] = 0
    indict['sub_unique'] = {}
    indict['sub_unique']['r1'] = 0
    indict['sub_unique']['r2'] = 0
    indict['sub_unique']['splits'] = 0
    indict['sub_unique']['split_types'] = {}

    for types in split_types:
        indict['sub_unique']['split_types'][types] = 0

    outdict = OrderedDict()
    outdict['templates'] = 0
    outdict['unmapped'] = 0
    outdict['multiple'] = 0
    outdict['sub_multiple'] = {}
    outdict['sub_multiple']['r1'] = 0
    outdict['sub_multiple']['r2'] = 0
    outdict['unique'] = 0
    outdict['sub_unique'] = {}
    outdict['sub_unique']['r1'] = 0
    outdict['sub_unique']['r2'] = 0
    outdict['sub_unique']['splits'] = 0

    metrics.append(indict)
    metrics.append(outdict)
    return metrics


def collect(sam_hits, samdict, ID):
    if ID is None:
        ID = samdict['qname']

    if ID == samdict['qname']:
        sam_hits.append(samdict)
        process_switch = 0
    else:
        process_switch = 1
        ID = samdict['qname']

    return sam_hits, samdict, ID, process_switch


def group_sam_hits_by_mappings(sam_hits):
    """
    Returns  a dictionary with:

        dict->mappingnumber->[0] = list
        dict->mappingnumber->[1] = list

        0 = first read, 1 = second read

    """

    # see if paired end or single and no mixture, multiple

    aln_dict = dict()

    for samdict in sam_hits:
        if 'XI' not in samdict['opt']:
            # pp.pprint(samdict)
            raise Exception('XI flag missing')

        XI_val = samdict['opt']['XI'][0]

        if XI_val not in aln_dict:
            for i in [0, 1]:
                aln_dict.setdefault(XI_val, {})[i] = []

        if args.read_type == 'single':
            samdict['flag'] = setBit(samdict['flag'], 6)
            samdict['flag'] = clearBit(samdict['flag'], 7)

        fstsnd = first_or_second_read(samdict['flag'])
        aln_dict[XI_val][fstsnd].append(samdict)

    return aln_dict


def process_split_reads(aln_dict):
    """ wrapper function for check_split_read and combine_split_read """
    for XI, template in aln_dict.items():
        for i in [0, 1]:
            list_of_fragments = template[i]
            tlength = len(list_of_fragments)
            if tlength > 1:
                info = check_split_read(list_of_fragments)
                template[i] = combine_split_reads(list_of_fragments, info)
    return aln_dict


def correct_flags_and_mate_information(aln_dict):
    """ horrible section needs to be cleaned up """

    suitable = []
    single_read_placement = dict()

    for XI, template in aln_dict.items():
        n_fst = len(template[0])
        n_snd = len(template[1])

        if args.read_type == 'single':
            if (n_fst == 1 and n_snd == 0):
                suitable.append(XI)
                single_read_placement[XI] = 0
            elif (n_fst == 0 and n_snd == 1):
                suitable.append(XI)
                single_read_placement[XI] = 1

        if args.read_type == 'paired':
            if (n_fst == 1 and n_snd == 1):
                suitable.append(XI)
            elif (n_fst == 1 and n_snd == 0):
                suitable.append(XI)
                single_read_placement[XI] = 0
            elif (n_fst == 0 and n_snd == 1):
                suitable.append(XI)
                single_read_placement[XI] = 1

    for i in suitable:
        template = aln_dict[i]
        NH = len(aln_dict.keys())

        if args.read_type == 'single':
            samdict = template[single_read_placement[i]][0]
            samdict = samdict_set_NH(samdict, NH)
            samdict = samdict_set_mapq(samdict)
            samdict = add_XS_field(samdict, [])

            aln_dict[i][single_read_placement[i]][0] = samdict

        if args.read_type == 'paired':
            # one read maps the other is unmapped
            if i in single_read_placement.keys():

                samdict = template[single_read_placement[i]][0]
                samdict = samdict_set_NH(samdict, NH)
                samdict = samdict_set_mapq(samdict)
                samdict = samdict_set_mate_unmapped(samdict)
                aln_dict[i][single_read_placement[i]][0] = samdict
                samdict = add_XS_field(samdict, [])

            else:

                samdicts = [template[0][0], template[1][0]]
                for n in [0, 1]:
                    samdicts[n] = samdict_set_NH(samdicts[n], NH)
                    samdicts[n] = samdict_set_mapq(samdicts[n])
                    samdicts[n] = add_XS_field(samdicts[n], [])

                (template[0][0], template[1][0]) = samdict_set_pe_mate(
                    samdicts[0], samdicts[1])
                aln_dict[i] = template

    aln_dict['suitable'] = suitable
    aln_dict['single_read_placement'] = single_read_placement

    return aln_dict


def collect_metrics(aln_dict, metrics):
    """
    Counts number of templates: example:
    This section needs a major overhaul, to many if and elses
    """

    indict = metrics[0]
    outdict = metrics[1]

    XI = aln_dict.keys()
    XI.remove('suitable')
    XI.remove('single_read_placement')
    XI_len = len(XI)

    suitable = aln_dict['suitable']

    if XI_len == 0:
        raise Exception(
            "Jesus Christ empty aln dict should not rise from the dead!")

    first = aln_dict[0][0]
    second = aln_dict[0][1]

    indict['templates'] += 1
    if XI_len > 1:
        indict['multiple'] += 1
        if len(first) > 0:
            indict['sub_multiple']['r1'] += 1
        if len(second) > 0:
            indict['sub_multiple']['r2'] += 1

        if len(suitable) > 0:
            outdict['templates'] += 1
            outdict['multiple'] += 1
            if len(first) > 0:
                outdict['sub_multiple']['r1'] += 1
            if len(second) > 0:
                outdict['sub_multiple']['r2'] += 1

    else:
        indict['unique'] += 1
        if len(suitable) > 0:
            outdict['unique'] += 1
            outdict['templates'] += 1

        if len(first) > 0:
            indict['sub_unique']['r1'] += 1
            split_type = is_split(first)

            if (split_type):
                indict['sub_unique']['splits'] += 1
                indict['sub_unique']['split_types'][split_type] += 1

            if len(suitable) > 0:
                outdict['sub_unique']['r1'] += 1

            if (split_type == 'OK'):
                outdict['sub_unique']['splits'] += 1

        if len(second) > 0:
            indict['sub_unique']['r2'] += 1

            split_type = is_split(second)
            if (split_type):
                indict['sub_unique']['splits'] += 1
                indict['sub_unique']['split_types'][split_type] += 1

            if len(suitable) > 0:
                outdict['sub_unique']['r2'] += 1

            if (split_type == 'OK'):
                outdict['sub_unique']['splits'] += 1

    return (aln_dict, [indict, outdict])


def output(aln_dict):
    single_read_placement = aln_dict['single_read_placement']
    suitable = aln_dict['suitable']

    if len(suitable) == 0:
        return

    for i in suitable:
        template = aln_dict[i]
        tlist = [0, 1]
        if i in single_read_placement.keys():
            tlist = [single_read_placement[i]]

        for read in tlist:
            sam_line = make_sam_line(template[read][0])
            args.outfile.write(sam_line + '\n')


def output_splits(aln_dict):
    XI = aln_dict.keys()
    XI.remove('suitable')
    XI.remove('single_read_placement')
    XI_len = len(XI)

    for entry in XI:
        template = aln_dict[entry]
        for i in [0, 1]:
            tlist = template[i]
            if len(tlist) > 1:
                for samdict in tlist:
                    sam_line = make_sam_line(samdict)
                    sys.stderr.write(sam_line)
                    args.outfile_splits.write(sam_line + '\n')


# called by functions used in main


def read_sam_line(line):
    """
    Returns the SAM line in a dictionary with optional fields
    stored in a seprate dictionary under samdict['opt']
    """
    fields = line.split('\t')
    keylist = 'qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual'.split(
        ', ')

    optional_fields = fields[11:]
    samdict = dict(zip(keylist, fields[0:11]))

    for numbers in ['flag', 'pos', 'mapq', 'isize']:
        samdict[numbers] = int(samdict[numbers])

    for opt_entry in optional_fields:
        (opt_id, opt_def, opt_val) = opt_entry.split(':')
        if opt_def == 'i':
            opt_val = int(opt_val)
        if opt_def == 'f':
            opt_val = float(opt_val)
        samdict.setdefault('opt', {})[opt_id] = [opt_val, opt_def]

    return samdict


def make_sam_line(samdict):
    """
    Takes a samdict and glues all the variables together to produce a SAM line
    """

    keylist = 'qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual'.split(
        ', ')
    pre_line = []

    # quick fix for single
    if args.read_type == 'single':
        samdict['flag'] = clearBit(samdict['flag'], 6)

    for key in keylist:
        pre_line.append(samdict[key])

    opt_dict = samdict['opt']
    for key in sorted(opt_dict.iterkeys()):
        tlist = [key,
                 opt_dict[key][1],
                 str(opt_dict[key][0])
                 ]
        opt = ':'.join(tlist)
        pre_line.append(opt)

    samline = "\t".join(str(v) for v in pre_line)

    return samline


def check_split_read(list_of_fragments):
    '''
    Takes a list of samdicts returns info dict with combinatorial status
    Insures that the fragment of split reads can be combined into one SAM line.
    The following checks are carried out:

        - Same chromosome
        - Same strand orientatien +/-
        - userdefined maximum distance allowed between fragments
        - No overlap between fragments
        - No switchorder within fragments
    '''

    info = dict()
    info['type'] = 'OK'

    # run checks
    # all frags should be on the same chr
    tlist = []
    for samdict in list_of_fragments:
        tlist.append(samdict['rname'])

    if not all_same(tlist):
        info['type'] = 'not_same_chr'

    # all frags should be on the same strand
    if info['type'] == 'OK':
        tlist = []
        for samdict in list_of_fragments:
            tlist.append(testBit(samdict['flag'], 4))

        if not all_same(tlist):
            info['type'] = 'not_same_strand'

    # all frags should be in args.distance
    # all frags should not overlap with each other
    list_of_fragments_by_pos = sorted(
        list_of_fragments, key=lambda k: k['pos'])
    if info['type'] == 'OK':
        tlist_start = []
        tlist_end = []
        for samdict in list_of_fragments_by_pos:
            tlist_start.append(samdict['pos'])
            seq_length = get_genomic_sequence_length(samdict['cigar'])
            tlist_end.append(samdict['pos'] + seq_length)

        tlist_start.pop(0)
        tlist_end.pop()

        for start, end in zip(tlist_start, tlist_end):
            diff = start - end
            if diff > args.distance:
                info['type'] = 'max_dist_breached'
                break
            if diff <= 0:
                info['type'] = 'fragment_overlap'

    # all frags should be in the right order
    if info['type'] == 'OK':

        tlist = []
        for samdict in list_of_fragments_by_pos:
            tlist.append(samdict['opt']['XQ'][0])

        samdict = list_of_fragments_by_pos[0]
        some_number = testBit(samdict['flag'], 4)

        if some_number > 0:
            tlist.reverse()

        if not is_sorted(tlist):
            info['type'] = 'order_switch'

    return info


def combinde_md_field(new_samdict, samdict):
    """
    Combines the MD field of two adjacent fragemnts.
    Segemehl seems to set always an integer at the start and end of the field.
    If no match is given 0 is expected.
    The fields are combined in partioning the first/previous  fragment (pf)
    into a  remainder  and the integer at the end of read for the next framgment (nf)
    the first intger is matched and the remainder of kept.
    """
    md_new = []

    md_pf = new_samdict['opt']['MD'][0]
    md_nf = samdict['opt']['MD'][0]

    prog = re.compile(r'\d+')
    pattern_pf = prog.finditer(md_pf)
    # md fields in segemehl seem to always start/end  with an integer, if a non match is encountered 0 is expected
    # there should be a number in there somewhere
    if pattern_pf is None:
        raise Exception(
            'No integer found in MD: {0} of previous Fragment \nreadname: {1}'.format(
                md_pf, new_samdict['qname']))

    # get last match through iterator
    last_match_pf = None
    for last_match_pf in pattern_pf:
        pass

    # is last match really at the end
    if last_match_pf.end() != len(md_pf):
        raise Exception(
            'No integer found at  end of MD: {0} of previous Fragment \nreadname: {1}'.format(
                md_pf, new_samdict['qname']))

    pattern_nf = prog.match(md_nf)
    if pattern_nf is None:
        raise Exception(
            'No integer found in MD: {0} of next Fragment \nreadname: {1}'.format(
                md_nf, samdict['qname']))

    if pattern_nf.start() > 0:
        raise Exception(
            'MD: {0} of next Fragment starts not with integer in \nreadname: {1}'.format(
                md_nf, samdict['qname']))

    # now go through the cases
    # keeping remains of md prev and next
    remainder_pf = md_pf[0:last_match_pf.start()]
    remainder_nf = md_nf[pattern_nf.end():]
    # for convienence
    pf = int(last_match_pf.group())
    nf = int(pattern_nf.group())

    #pp.pprint(['pf',md_pf, remainder_pf, pf])
    #pp.pprint(['nf',md_nf, remainder_nf, nf])

    # case both have matches
    if (pf > 0 and nf > 0):
        tmatch = str(pf + nf)
        md_new = ''.join([remainder_pf, tmatch, remainder_nf])
    # pf has matches nf something else , discard 0 of nf
    elif(pf >= 1 and nf == 0):
        md_new = ''.join([remainder_pf, str(pf), remainder_nf])
    # pf has matches nf something else , discard 0 of nf
    elif(pf == 0 and nf >= 1):
        md_new = ''.join([remainder_pf, str(nf), remainder_nf])
    # both have no matches so just discarding 0 of pf should be enough
    elif(pf == 0 and nf == 0):
        md_new = ''.join([remainder_pf, str(nf), remainder_nf])
    else:
        raise Exception(
            'pf {2} nf {3} unexpected MD case:\n{0}\n {1}'.format(
                new_samdict, samdict, pf, nf))

    return md_new


def combine_split_reads(list_of_fragments, info):
    """
    Combines the split read if possible in any case my YS: flag is appended with the state of the read
    The first fragment serves as template for the combined read in which
    """
    if info['type'] != 'OK':
        # mark not combinable split read in samdict
        for samdict in list_of_fragments:
            samdict.setdefault('opt', {})['YS'] = [info['type'], 'Z']
        return list_of_fragments
    else:
        list_of_fragments_by_pos = sorted(
            list_of_fragments, key=lambda k: k['pos'])
        new_samdict = list_of_fragments_by_pos[0]

        # need to set
        # nucleotide sequence ['seq']
        # mapping quality ['mapq'], dont know if all the frags have the same quality I just take the mean
        # base qualities['qual']
        ##MD  ['opt']['MD'][0]
        ##NH  ['opt']['NM'][0]

        mapq = []
        mapq.append(new_samdict['mapq'])
        XS_list = []
        for i, samdict in enumerate(list_of_fragments_by_pos[1:]):
         #       pp.pprint(i)
          #      pp.pprint(samdict)

            new_samdict['seq'] += samdict['seq']
            new_samdict['qual'] += samdict['qual']
            new_samdict['opt']['NM'][0] += samdict['opt']['NM'][0]

            seq_length = get_genomic_sequence_length(new_samdict['cigar'])
            diff = samdict['pos'] - (new_samdict['pos'] + seq_length)
            new_samdict['cigar'] = ''.join(
                [new_samdict['cigar'], str(diff), 'N', samdict['cigar']])

            if (args.library_type == 'fr-unstranded' and args.seq_type == 'RNA'):
                record_dict = args.genome_dict
                rname = new_samdict['rname']
                spliceA_start = (new_samdict['pos'] + seq_length - 1)
                spliceA_end = spliceA_start + 2
                spliceB_start = samdict['pos'] - 3
                spliceB_end = samdict['pos'] - 1

                subseq = str(record_dict[rname].seq[spliceA_start:spliceA_end])
                subseq += str(record_dict[rname].seq[spliceB_start:spliceB_end])

                # splice sites are found on the plus strand so transcript comes
                # from +
                if(subseq == "GTAG"or subseq == "GCAG" or subseq == "ATAC"):
                    XS_list.append('+')

                    # splice sites are found on the minus strand so transcript
                    # comes from -
                elif(subseq == "CTAC"or subseq == "CTGC" or subseq == "GTAT"):
                    XS_list.append('-')
                    # not a common splice site used so making no inferences
                else:
                    XS_list.append('dunno')

            mapq.append(samdict['mapq'])
            # MD fields suck  ....
            new_samdict['opt']['MD'][0] = combinde_md_field(
                new_samdict, samdict)

        new_samdict = clean_dict(new_samdict)

        new_samdict = add_XS_field(new_samdict, XS_list)

        mapq_mean = numpy.mean(mapq)
        new_samdict['mapq'] = int(numpy.floor(mapq_mean))
        new_samdict.setdefault('opt', {})['YS'] = [info['type'], 'Z']
        # print('samdict')

        return [new_samdict]


def add_XS_field(samdict, XS_list):

    # 0 = fst 1 =snd
    fstsnd = str(first_or_second_read(samdict['flag']))
    # 0 = sense 1 = reverse
    direction = str(read_orientation(samdict['flag']))

    XS_type = ''.join([fstsnd, direction])
    XS_dict = dict()
    XS_dict['00'] = '-'
    XS_dict['01'] = '+'
    XS_dict['10'] = '+'
    XS_dict['11'] = '-'
    # there  are sure cleaner ways

    XS_dict['+'] = '-'
    XS_dict['-'] = '+'

    XS_val = XS_dict[XS_type]

    if (args.library_type == 'fr-unstranded' and args.seq_type == 'RNA'):
        if (len(XS_list) == 0):
            pass
        elif ('dunno' in XS_list):
            pass
        elif not all_same(XS_list):
            pass
        else:
            samdict.setdefault('opt', {})['XS'] = [XS_list[0], 'A']

    elif (args.library_type == 'fr-firststrand' and args.seq_type == 'RNA'):
        # by definition
        # first read reverse -> comes from +
        # first read not reverse -> comes from -
        # second read reverse -> comes from -
        # second not read reverse -> comes from +

        samdict.setdefault('opt', {})['XS'] = [XS_val, 'A']

    elif (args.library_type == 'fr-secondstrand' and args.seq_type == 'RNA'):
        # by definition
        # first read reverse -> comes from -
        # first read not reverse -> comes from +
        # second read reverse -> comes from +
        # second not read reverse -> comes from -
        # since its the opposite outcome of fr-firstrand +/- are swapped

        samdict.setdefault('opt', {})['XS'] = [XS_dict[XS_val], 'A']

    else:
        pass
    return samdict


def clean_dict(samdict):
    """
    Removes unecessary private TAGS mostly split read specific.
    Oh TAG XA not found in manual where do you come from what do you mean
    for now I will swipe you clean
    """

    tags = ['XA', 'XX', 'XY', 'XQ', 'XL', 'XP', 'XU', 'XS', 'XC', 'XV', 'XT']
    for tag in tags:
        if tag in samdict['opt']:
            del samdict['opt'][tag]

    return samdict


def samdict_set_mate_unmapped(samdict):
    samdict['flag'] = setBit(samdict['flag'], 3)
    samdict['mrnm'] = '*'
    samdict['mpos'] = 0
    samdict['isize'] = 0
    return samdict


def samdict_set_pe_mate(samdict_A, samdict_B):
    """
    Sets:

        - mate reference name
        - mate pos
        - insertion size
        - proper pair (same chr && reverse read behind not reverse read)


    """

    tdict = dict()

    if (samdict_A['rname'] == samdict_B['rname']):
        samdict_A['mrnm'] = '='
        samdict_B['mrnm'] = '='

        # set isize
        if samdict_A['pos'] < samdict_B['pos']:
            tdict['left'] = samdict_A
            tdict['right'] = samdict_B
        else:
            #tdict.setdefault('right', default={})
            tdict['right'] = samdict_A
            tdict['left'] = samdict_B

        seq_length = get_genomic_sequence_length(tdict['right']['cigar'])
        isize = tdict['right']['isize'] + seq_length - tdict['left']['isize']

        tdict['right']['isize'] = isize
        tdict['left']['isize'] = isize * (-1)

    else:
        samdict_A['mrnm'] = samdict_B['rname']
        samdict_B['mrnm'] = samdict_A['rname']

        samdict_A['isize'] = 0
        samdict_B['isize'] = 0

    samdict_A['mpos'] = samdict_B['pos']
    samdict_B['mpos'] = samdict_A['pos']

    # if mate reverse set mate reverse flag in my info
    if testBit(samdict_B['flag'], 4) > 0:
        samdict_A['flag'] = setBit(samdict_A['flag'], 5)

    if testBit(samdict_A['flag'], 4) > 0:
        samdict_B['flag'] = setBit(samdict_B['flag'], 5)

    # proper pair
    # same chr read orientation must face each other ---> <---
    samdict_A['flag'] = clearBit(samdict_A['flag'], 1)
    samdict_B['flag'] = clearBit(samdict_B['flag'], 1)

    # same chr:
    if (samdict_A['rname'] == samdict_B['rname']):
        one_reverse = None
        # reverse read sholud always be after unreversed read
        if (testBit(samdict_A['flag'], 4) >
                0 and testBit(samdict_A['flag'], 5) == 0):
            if (samdict_B['pos'] <= samdict_A['pos']):
                one_reverse = 1

        elif (testBit(samdict_A['flag'], 4) == 0 and testBit(samdict_A['flag'], 5) > 0):
            if (samdict_A['pos'] <= samdict_B['pos']):
                one_reverse = 1

        if (one_reverse):
            samdict_A['flag'] = setBit(samdict_A['flag'], 1)
            samdict_B['flag'] = setBit(samdict_B['flag'], 1)

    return (samdict_A, samdict_B)


def get_genomic_sequence_length(cigar):
    """
    Returns sequence length in genomic space by
    doing the cigar opeartions.
    Does not catch malformed cigarlines
    """
    cigar_strings = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
    cigar_dict = dict.fromkeys(cigar_strings, 0)

    tlist = re.findall(r'\d+\D', cigar, flags=0)
    for i in tlist:
        cigar_type = i[-1:]
        count = int(i[:-1])
        cigar_dict[cigar_type] += count

    list_sum_of_length = ['M', 'I', 'S', '=', 'X']
    seq_length = 0
    for i in list_sum_of_length:
        seq_length += cigar_dict[i]

    return seq_length


def first_or_second_read(flag):
    """returns 0 if its first sequenced read and 1 for the second sequenced read """

    fst = testBit(flag, 6)
    snd = testBit(flag, 7)

    if (fst > 0 and snd == 0):
        return 0
    elif(fst == 0 and snd > 0):
        return 1
    else:
        raise Exception(
            'Not unique flag assignment fst: {0} snd: {1}'.format(
                fst, snd))


def read_orientation(flag):
    """returns 0 if sense/+ and 1 for reversed """

    orient = testBit(flag, 4)
    if (orient > 0):
        return 1
    else:
        return 0


def samdict_set_mapq(samdict):
    NH = samdict['NH']
    if (NH == 1):
        samdict['mapq'] = 50
    elif (NH == 2):
        samdict['mapq'] = 5
    else:
        samdict['mapq'] = 0

    return samdict


def samdict_set_NH(samdict, NH):

    samdict['flag'] = clearBit(samdict['flag'], 8)
    samdict['NH'] = NH

    if (NH > 1):
        samdict['flag'] = setBit(samdict['flag'], 8)

    return samdict


def is_sorted(lst, key=lambda x, y: x < y):
    for i, el in enumerate(lst[1:]):
        if key(el, lst[i]):
            return False
    return True


def test_set_flag(flag, bit_to_check, response, bit_to_set):
    if testBit(flag, bit_to_check) > response:
        flag = setBit(flag, bit_to_set)
    return flag


def all_same(items):
    """
    tests if all items in alist are the same
    http://stackoverflow.com/a/3787983
    """
    return all(x == items[0] for x in items)


def testBit(int_type, offset):
    """    testBit() returns a nonzero result, 2**offset, if the bit at 'offset' is one.    """
    mask = 1 << offset
    return(int_type & mask)


def setBit(int_type, offset):
    """setBit() returns an integer with the bit at 'offset' set to 1."""
    mask = 1 << offset
    return(int_type | mask)


def clearBit(int_type, offset):
    """clearBit() returns an integer with the bit at 'offset' cleared."""
    mask = ~(1 << offset)
    return(int_type & mask)


def toggleBit(int_type, offset):
    """toggleBit() returns an integer with the bit at 'offset' inverted, 0 -> 1 and 1 -> 0"""
    mask = 1 << offset
    return(int_type ^ mask)


def is_split(samdicts):
    samdict = samdicts[0]

    if 'YS' in samdict['opt']:
        return samdict['opt']['YS'][0]
    else:
        return None


def output_metrics(metrics):
    args.logfile.write(yaml.dump(metrics, default_flow_style=None))


def check_sort_order_sam_file(lst):
    """NOT IMPLEMENTED samtools sort and pyhton sort do not get along at the moment"""
    pass
    # if sorted(lst) == lst:
    #    pass
    # else:
    #    pp.pprint(lst)
    #    pp.pprint(sorted(lst))
    #    raise StandardError('SAM file not sorted by name:\n\t fst: {0}\n\t snd: {1} entry'.format(lst[0], lst[1]))


def filter_for_snp_calling(aln_dict):

    XI = aln_dict.keys()
    XI.remove('suitable')
    XI.remove('single_read_placement')
    NH = len(XI)
    if NH > 1:
        aln_dict['suitable'] = []
        return aln_dict

    fst = aln_dict[0][0]
    snd = aln_dict[0][1]

    n_fst = len(fst)
    n_snd = len(snd)

    if n_fst == n_snd:
        if fst[0]['mrnm'] != '=':
            aln_dict['suitable'] = []
            return aln_dict

    return aln_dict


def main(args):
    ID = None
    # print(args)
    out = args.outfile
    sam_hits = []
    metrics = init_metrics()

    for line in args.infile:
        # in case the SAM header is piped
        if line[0] == '@':
            args.outfile.write(line)
            continue
        line = line.rstrip('\n')
        samdict = read_sam_line(line)

        # take all reads same name and put them into a list "sam_hits"
        (sam_hits, samdict, ID, process_switch) = collect(sam_hits, samdict, ID)

        # new readname  encountered start processesing
        if process_switch == 1:
            # check_sort_order_sam_file([sam_hits[0]['qname'],ID])
            # look at XI paires and determine if only paired or single are in
            # the collection
            aln_dict = group_sam_hits_by_mappings(sam_hits)
            aln_dict = process_split_reads(aln_dict)
            aln_dict = correct_flags_and_mate_information(aln_dict)

            if (args.filter_snp_calling):
                aln_dict = filter_for_snp_calling(aln_dict)
            aln_dict, metrics = collect_metrics(aln_dict, metrics)
            # pp.pprint((aln_dict))

            output(aln_dict)

            if (args.outfile_splits):
                output_splits(aln_dict)
            aln_dict = []
            sam_hits = []
            sam_hits.append(samdict)

    # last block after while iteration
    aln_dict = group_sam_hits_by_mappings(sam_hits)
    aln_dict = process_split_reads(aln_dict)
    aln_dict = correct_flags_and_mate_information(aln_dict)
    if (args.filter_snp_calling):
        aln_dict = filter_for_snp_calling(aln_dict)

    aln_dict, metrics = collect_metrics(aln_dict, metrics)
    output(aln_dict)
    if (args.outfile_splits):
        output_splits(aln_dict)
    output_metrics(metrics)


if __name__ == '__main__':
    args = read_arguments()
    eval_arguments(args)
    main(args)

__version__ = '0.001'
