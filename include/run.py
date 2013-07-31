import abstract_step
import misc
import os

class Run(object):
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
        if key in self._private_info and value != self._private_info[key]:
            raise StandardError(
                "You're trying to overwrite private info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._private_info[key]))
        self._private_info[key] = value

    def add_public_info(self, key, value):
        if key in self._public_info and value != self._public_info[key]:
            raise StandardError(
                "You're trying to overwrite public info %s with %s, "
                "but there's already a different value stored: %s." %
                (key, value, self._public_info[key]))
        self._public_info[key] = value

    def add_output_file(self, tag, out_path, in_paths):
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

    def output_files(self):
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
        temp = self.output_files()
        if len(temp[annotation]) != 1:
            raise StandardError("More than one output file declared for out/%s." % annotation)
        return temp[annotation].keys()[0]

    def get_output_files_for_annotation_and_tags(self, annotation, tags):
        temp = self.output_files()
        return misc.assign_strings(temp[annotation].keys(), tags)

    def get_input_files_for_output_file(self, out_path):
        temp = self.output_files()
        for tag in temp.keys():
            if out_path in temp[tag].keys():
                return sorted(temp[tag][out_path])
        raise StandardError("Sorry, your output file couldn't be found in the dictionary: %s." % out_path)

    def public_info(self, key):
        return self._public_info[key]

    def has_public_info(self, key):
        return (key in self._public_info)

    def private_info(self, key):
        return self._private_info[key]

    def as_dict(self):
        result = dict()
        result['output_files'] = self._output_files
        result['private_info'] = self._private_info
        result['public_info'] = self._public_info
        result['run_id'] = self._run_id
        return result
