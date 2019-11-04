import os
seq_pipeline_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/../python_env/bin/activate_this.py' % seq_pipeline_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import argparse

'''

'''

######################
### help functions ###
######################

def get_cleaned_line(line):
    line = line.split('\t')
    gene = line[0]
    count = line[-1].rstrip('\n')
    return [gene, float(count)]

def is_headline(line):
    if line[0] == '#':
        return True
    if line[0:6] == 'Geneid':
        return True
    return False

#####################
### featureCounts ###
#####################
def fc_merge(files):
    basic_data = []
    for i, infile in enumerate(files):
        # take the first as template
        if i == 0:
            with open(infile, 'r') as f:
                for i, line in enumerate(f):
                    # first and second line in featureCounts output file are headlines
                    if i == 0:
                        new_line = '# merged featureCounts output files: %s' % (
                        ' '.join(files))
                        basic_data.append(new_line)
                        continue

                    if i == 1:
                        # remove filename at the end
                        new_line = line.rstrip('\n').split('\t')
                        new_line = new_line[0:-1]
                        basic_data.append('\t'.join(new_line))
                        continue

                    basic_data.append(line.rstrip('\n').split('\t'))
            continue

        with open(infile, 'r') as f:
            for i, line in enumerate(f):
                # continue if headline
                if is_headline(line):
                    continue

                line = line.strip('\n').split('\t')
                basic_data[i][-1] = str(int(basic_data[i][-1]) + int(line[-1]))

    return basic_data

def fc_write_to_file(out_file, data):
    with open(out_file, 'w') as f:
        for i, line in enumerate(data):
            if i == 0 or i == 1:
                f.write(line + '\n')
                continue

            f.write('\t'.join(map(str, line)) + '\n')

###################
### htseq_count ###
###################
def htc_merge(files):
    basic_data = []
    for i, file in enumerate(files):
        if i == 0:
            with open(file, 'r') as f:
                for line in f:
                    basic_data.append(get_cleaned_line(line))
            continue

        with open(file, 'r') as f:
            for i, line in enumerate(f):
                line = get_cleaned_line(line)
                basic_data[i][1] += line[1]

    return basic_data


def htc_write_to_file(out_file, data):
    with open(out_file, 'w') as f:
        for i, line in enumerate(data):
            f.write('\t'.join(map(str, line)) + '\n')

###############
### merging ###
###############

def read_args():
    parser = argparse.ArgumentParser(
        description='merges genecounts from different tools (htseq_count, featureCounts)'
    )
    parser.add_argument(
        '--tool_name',
        '-t',
        nargs=1,
        help='tool name (htseq_count: htc, featureCounts: fc)'
    )
    parser.add_argument(
        '--outpath',
        '-o',
        nargs=1,
        help='output path'
    )
    parser.add_argument(
        '--out_file_pattern',
        '-p',
        nargs='?',
        help='file name pattern for output files, a counter will be added automatically at the end (<pattern>_1.fastq)'
    )

    # TODO: catch missing args
    args, other_args = parser.parse_known_args()

    return [args, other_args]

def main(args):
    options, raw_file_list = args
    file_list = raw_file_list[0].split(' ')

    tool_name = options.tool_name[0]
    outpath = options.outpath[0]
    if outpath != '' and outpath[-1] != '/':
        outpath += '/'

    out_file_pattern = options.out_file_pattern

    data = []
    if tool_name == 'fc':
        data = fc_merge(file_list)
        default_filename = 'featureCounts_merged.txt'
        # write data to file
        out_file = outpath + out_file_pattern if out_file_pattern is not None else outpath + default_filename
        fc_write_to_file(out_file, data)
    elif tool_name == 'htc':
        data = htc_merge(file_list)
        default_filename = 'htseq_count_merged.txt'
        # write data to file
        out_file = outpath + out_file_pattern if out_file_pattern is not None else outpath + default_filename
        htc_write_to_file(out_file, data)
    else:
        # TODO: throw exception
        pass


if __name__ == '__main__':
    args = read_args()
    main(args)
