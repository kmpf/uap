import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import argparse

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
        nargs='*',
        help='output file'
    )

    args = parser.parse_args()
    return args

def read_mapping_file(mapping_file):
    mapping_data = dict()
    with open(mapping_file, 'r') as f:
        for line in f:
            line_data = line.rstrip('\n').split('\t')
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

def write_output_file(output_file, result_data):
    with open(output_file, 'w') as f:
        for g_id, count in result_data.iteritems():
            f.write('%s\t%f\n' % (g_id, count))

def main(args):
    # TODO: mapping file exists?
    mapping_file = args.mapping_file[0]
    # TODO: input file exists? Check format?
    input_file = args.infile[0]
    # TODO: tool_name kallisto or salmon?
    tool_name = args.tool_name[0]
    output_file = args.outfile[0] if args.outfile else None

    mapping_data = read_mapping_file(mapping_file)
    result_data = read_input_file(input_file, mapping_data, tool_name)
    write_output_file(output_file, result_data)

if __name__ == '__main__':
    args = read_args()
    main(args)