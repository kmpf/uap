import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import argparse
from collections import defaultdict
import sys
import json
'''
calculates gene counts from transcript counts
tcount2gcounts.py -m <transcript_gene_mapping.csv> -i <transcript_counts.csv> -o <gene_counts.csv>
'''

def read_args():
    parser = argparse.ArgumentParser(
        description='calculates gene counts from transcript counts'
    )
    parser.add_argument(
        '--mapping-file',
        '-m',
        type=str,
        nargs=1,
        help='transcript to gene mapping file. Required Format example (per row): ' \
             'ENST00000527779.1\tENSG00000137692.11'
    )
    parser.add_argument(
        '--infile',
        '-i',
        nargs=1,
        help='transcript counts input file.'
    )
    parser.add_argument(
        '--tool_name',
        '-t',
        nargs=1,
        help='source tool of input file. Possible values: kallisto, salmon'
    )
    parser.add_argument(
        '--outfile',
        '-o',
        type=argparse.FileType('w'),
        nargs='?',
        help='output file',
        default=sys.stdout
    )


    parser.add_argument(
        '--kallisto-extended',
        action='store_true',
        help='summarizes tpm and counts sets length to NA'

    )

    args = parser.parse_args()
    return args

def read_mapping_file(mapping_file):
    mapping_data = dict()
    
    with open(mapping_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue 

            line_data = line.rstrip('\n').split('\t')
            if mapping_file.endswith('gtf'):
                data = line_data[8].split()
                if all(elem in data for elem in ['gene_id', 'transcript_id']):
                    d = dict(zip(data[::2], data[1::2]))
                    t_id = d['transcript_id'].rstrip(';').replace('"', '')
                    g_id = d['gene_id'].rstrip(';').replace('"', '')
                else:
                    continue
            
                    
            else:
                t_id = line_data[0]
                g_id = line_data[1]
            mapping_data[t_id] = g_id

    # TODO: mapping data empty after extracting?
    return mapping_data

def read_input_file(input_file, mapping_data, tool_name):
    result_data = dict()

    count_index = 3 if tool_name == 'kallisto' else 4

    with open(input_file, 'r') as f:
        for lnum, line in enumerate(f):
            if lnum != 0:
                line_data = line.rstrip('\n').split('\t')
                transcript = line_data[0].strip('"')
                count = line_data[count_index]

                if mapping_data[transcript] not in result_data:
                    result_data[mapping_data[transcript]] = 0

                result_data[mapping_data[transcript]] += float(count)

    return result_data


def read_input_file_k(input_file, mapping_data):
    res =  defaultdict(dict)

    info = ['t_target_id', 't_length', 't_eff_length', 't_est_count', 't_tpm'] 

    with open(input_file, 'r') as f:
        for lnum, line in enumerate(f):
            
            if lnum != 0:
                data  = line.strip('"').rstrip('\n').split('\t')

                transcript = data[0]

                if mapping_data[transcript] not in res:
                    for i in range(5):
                        res[mapping_data[transcript]][info[i]] = str(data[i])
                    res[mapping_data[transcript]]['target_id'] = str(mapping_data[transcript])
                    res[mapping_data[transcript]]['length'] = 'NA'
                    res[mapping_data[transcript]]['eff_length'] = 'NA'
                    res[mapping_data[transcript]]['est_count'] = float(data[3])
                    res[mapping_data[transcript]]['tpm'] = float(data[4])
                        
                else:
                    for i in range(5):
                        res[mapping_data[transcript]][info[i]] += ',' + str(data[i])
                    res[mapping_data[transcript]]['est_count'] += float(data[3])
                    res[mapping_data[transcript]]['tpm'] += float(data[4])


    return res



def write_output_file(output_file, result_data):
    for g_id, count in result_data.items():
        output_file.write('%s\t%f\n' % (g_id, count))

def write_output_file_k(output_file, result_data):

    header = ['target_id', 'length', 'eff_length', 'est_count', 'tpm',
              't_target_id', 't_length', 't_eff_length', 't_est_count', 't_tpm'] 
    
    output_file.write("\t".join(header) + '\n')
    for g_id, entry in result_data.items():
        out = [str(entry[x]) for x in header]
        output_file.write("\t".join(out) + "\n")



def main(args):
    # TODO: mapping file exists?
    mapping_file = args.mapping_file[0]
    # TODO: input file exists? Check format?
    input_file = args.infile[0]
    # TODO: tool_name kallisto or salmon?
    tool_name = args.tool_name[0]


    mapping_data = read_mapping_file(mapping_file)

    if args.kallisto_extended:
        result_data = read_input_file_k(input_file, mapping_data)
        write_output_file_k(args.outfile, result_data)
    else:
        result_data = read_input_file(input_file, mapping_data, tool_name)
        write_output_file(args.outfile, result_data)

if __name__ == '__main__':
    args = read_args()
    main(args)
