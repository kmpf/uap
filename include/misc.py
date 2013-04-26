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
