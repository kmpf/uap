import hashlib
import json
import os
import re

def assign_strings(paths, tags):
    '''
    Assign strings (path names, for example) to tags. Example:
    
    - paths = ['RIB0000794-cutadapt-R1.fastq.gz', 'RIB0000794-cutadapt-R2.fastq.gz']
    - tags = ['R1', 'R2']
    - result = { 'R1': 'RIB0000794-cutadapt-R1.fastq.gz', 'R2': 'RIB0000794-cutadapt-R2.fastq.gz' }
      
    If this is not possible without ambiguities, a StandardError is thrown.
    '''

    def check_candidate(paths, tags, head, tail):
        chopped = []
        for path in paths:
            if path[:len(head)] != head:
                return None
            if path[-len(tail):] != tail:
                return None
            chopped.append((path[len(head):-len(tail)], path))

        if [_[0] for _ in sorted(chopped)] == sorted(tags):
            result = {}
            for _ in sorted(chopped):
                result[_[0]] = _[1]
            return result

        return None


    results = {}
    if len(paths) != len(tags):
        raise StandardError("Number of tags must be equal to number of paths")
    
    for tag in tags:
        for path in paths:
            result_candidate = {}
            if tag in path:
                # find all occurences of tag in path
                offset = 0
                while path.find(tag, offset) >= 0:
                    index = path.find(tag, offset)
                    head = path[:index]
                    tail = path[(index+len(tag)):]
                    # now try chopping off head and tail from every path
                    # and see whether we can unambiguously assign a path
                    # to every tag, if yes, we have a result candidate
                    result_candidate = check_candidate(paths, tags, head, tail)
                    if result_candidate:
                        results[json.dumps(result_candidate, sort_keys = True)] = result_candidate
                    offset = index + 1
                    
    if len(results) != 1:
        raise StandardError("Unable to find an unambiguous mapping.")
    
    return results[results.keys()[0]]

def natsorted(l):
    '''
    Return a 'naturally sorted' permutation of l.
    
    Credits: http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html
    '''
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key=alphanum_key)
    
def str_to_sha1(s):
    return hashlib.sha1(s).hexdigest()

def bytes_to_str(num):
    for _, x in enumerate(['bytes','k','M','G']):
        if num < 1024.0:
            if _ == 0:
                return "%d %s" % (num, x)
            else:
                return "%1.1f %sB" % (num, x)
        num /= 1024.0
    return "%1.1f %sB" % (num, 'T')

def duration_to_str(duration, long = False):
    value = str(duration)
    if not long:
        if 'days' in value:
            value = value.replace(' days,', 'd')
        if 'day' in value:
            value = value.replace(' day,', 'd')
        if 'd' in value and ':' in value and (value.index(':') - value.index('d')) != 4:
            value = value[:value.index('d') + 1] + ' ' + value[value.index('d') + 1:]
    if '.' in value:
        value = value[0:value.index('.') + 2]
    return value

def append_suffix_to_path(path, suffix):
    dirname, filename = os.path.split(path)
    if '.' in filename:
        basename = filename[:filename.index('.')]
        extension = filename[filename.index('.'):]
    else:
        basename = filename
        extension = ''
    filename = basename + '-' + suffix + extension
    return os.path.join(dirname, filename)