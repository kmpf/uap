import abstract_step
import misc
import os

class Run(object):
    '''
    The Run class is a helper class which represents a run in a step. Declare runs
    inside AbstractStep.declare_runs() via::
    
        with self.declare_run(run_id) as run:
            # declare output files, private and public info here
            
    After that, use the available methods to configure the tun.
    '''
    def __init__(self, step, run_id):
        if '/' in run_id:
            raise StandardError("Error: A run ID must not contain a slash: %s." % run_id)
        self._step = step
        self._run_id = run_id
        self._private_info = dict()
        self._public_info = dict()
        self._output_files = dict()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def add_private_info(self, key, value):
        '''
        Add private information to a run. Use this to store data which you will need 
        when the run is executed. As opposed to public information, private information 
        is not visible to subsequent steps.
        
        You can store paths to input files here, but not paths to output files as
        their expected location is not defined until we're in *AbstractStep.execute*
        (hint: they get written to a temporary directory inside *execute()*).
        '''
        if key in self._private_info and value != self._private_info[key]:
            raise StandardError(
                "You're trying to overwrite private info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._private_info[key]))
        self._private_info[key] = value

    def add_public_info(self, key, value):
        '''
        Add public information to a run. For example, a FASTQ reader may store the
        index barcode here for subsequent steps to query via 
        ``AbstractStep.find_upstream_info()``.
        '''
        if key in self._public_info and value != self._public_info[key]:
            raise StandardError(
                "You're trying to overwrite public info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._public_info[key]))
        self._public_info[key] = value

    def add_output_file(self, tag, out_path, in_paths):
        '''
        Add an output file to this run. Output file names must be unique across all
        runs defined by a step, so it may be a good idea to include the run_id into
        the output filename.

        - *tag*: You must specify the connection annotation which must have been previously 
          declared via *AbstractStep.add_connection("out/...")*, but this doesn't 
          have to be done in the step constructor, it's also possible in *declare_runs()*
          right before this method is called.
        
        - *out_path*: The output file path, without a directory. The pipeline assigns
          directories for you (this parameter must not contain a slash).
        
        - *in_paths*: A list of input files this output file depends on. It is 
          **crucial** to get this right, so that the pipeline can determine which steps
          are up-to-date at any given time. You have to specify absolute paths here,
          including a directory, and you can obtain them via 
          *AbstractStep.run_ids_and_input_files_for_connection* and related functions.
        '''
        
        # make sure there's no slash in out_path unless it's a source step
        if '/' in out_path and abstract_step.AbstractSourceStep not in self._step.__class__.__bases__:
            raise StandardError("There must be no slash (/) in any output "
                "file declared by a step: %s." % out_path)
        # make sure tag was declared with an outgoing connection
        if 'out/' + tag not in self._step._connections:
            raise StandardError("Invalid output_file tag '%s' in %s. "
                "You might want to add self.add_connection('out/%s') "
                "to the constructor of %s."
                % (tag, str(self._step), tag, self._step.__module__))

        if tag not in self._output_files:
            self._output_files[tag] = dict()

        if out_path in self._output_files[tag]:
            raise StandardError(
                "You're trying to re-add an output file which has already "
                "been declared: %s." % out_path)

        self._output_files[tag][out_path] = in_paths


    def get_output_files(self):
        '''
        Return a dictionary of all defined output files, grouped by connection 
        annotation::
        
           annotation_1:
               out_path_1: [in_path_1, in_path_2, ...]
               out_path_2: ...
           annotation_2: ...
        '''
        result = dict()
        for tag in self._output_files:
            result[tag] = dict()
            for out_path, in_paths in self._output_files[tag].items():
                directory = self._step.get_output_directory_du_jour()
                if directory != None:
                    full_path = os.path.join(directory, out_path)
                else:
                    full_path = out_path
                result[tag][full_path] = in_paths
        return result

    def get_single_output_file_for_annotation(self, annotation):
        '''
        Retrieve exactly one output file of the given annotation, and crash
        if there isn't exactly one.
        '''
        temp = self.get_output_files()
        if len(temp[annotation]) != 1:
            raise StandardError("More than one output file declared for out/%s." % annotation)
        return temp[annotation].keys()[0]

    def get_output_files_for_annotation_and_tags(self, annotation, tags):
        '''
        Retrieve a set of output files of the given annotation, assigned to
        the same number of specified tags. If you have two 'alignment' output files
        and they are called *out-a.txt* and *out-b.txt*, you can use this function like this:
        
        - *tags*: ['a', 'b']
        - result: {'a': 'out-a.txt', 'b': 'out-b.txt'}
        '''
        temp = self.get_output_files()
        return misc.assign_strings(temp[annotation].keys(), tags)

    def get_input_files_for_output_file(self, out_path):
        '''
        Return all input files a given output file depends on.
        '''
        temp = self.get_output_files()
        for tag in temp.keys():
            if out_path in temp[tag].keys():
                return sorted(temp[tag][out_path])
        raise StandardError("Sorry, your output file couldn't be found in the dictionary: %s." % out_path)

    def get_public_info(self, key):
        '''
        Query public information which must have been previously stored via *add_public_info()*.
        '''
        return self._public_info[key]

    def has_public_info(self, key):
        '''
        Query whether a piece of public information has been defined.
        '''
        return (key in self._public_info)

    def get_private_info(self, key):
        '''
        Query private information which must have been previously stored via *add_private_info()*.
        '''
        return self._private_info[key]

    def has_private_info(self, key):
        '''
        Query whether a piece of public information has been defined.
        '''
        return (key in self._private_info)

    def as_dict(self):
        result = dict()
        result['output_files'] = self._output_files
        result['private_info'] = self._private_info
        result['public_info'] = self._public_info
        result['run_id'] = self._run_id
        return result
