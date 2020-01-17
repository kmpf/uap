#!/usr/bin/env python

import logging
import glob
import os
import yaml
import time
import pipeline

'''
This script
'''

logger = logging.getLogger("uap_logger")

def write_output_file(data, file_path):
    f = open(file_path, 'w')

    for line in data:
        f.write(','.join(map(str, line)) + '\n')

    f.close()

def main(args):
    p = pipeline.Pipeline(arguments=args)

    # Compile a list of all tasks
    task_list = list()
    # Only use tasks listed in args.run
    arguments = vars(args)

    if arguments['run_id'] and len(arguments['run_id']) >= 1:
        for task_id in arguments['run_id']:
            if '/' in task_id:
                task = p.task_for_task_id[task_id]
                task_list.append(task)
            else:
                for task in p.all_tasks_topologically_sorted:
                    if str(task)[0:len(task_id)] == task_id:
                        task_list.append(task)
    # or take all available tasks
    else:
        task_list = p.all_tasks_topologically_sorted


    output_data = [
        ['# uap runtime-info', time.strftime("%c"), p.git_hash_tag],
        ['# uap  config', p.get_config_filepath()],
        ['# step', 'CPUs [%]', 'requested CPUs', 'RAM [MB]', 'Duration [s]', 'total start time', 'total end time'],
        ['step', 'cpu_percent', 'requested_cpu', 'ram_mb', 'duration_seconds', 'total_start_time', 'total_end_time']
    ]


    for task in task_list:
        outdir = task.get_run().get_output_directory()
        anno_files = glob.glob(os.path.join(
            outdir, ".%s*.annotation.yaml" % task.get_run().get_run_id()
        ))

        yaml_files = {os.path.realpath(f) for f in anno_files \
                      if os.path.islink(f)}

        for y in yaml_files:
            annotation_data = Annotation_Data(y)

            output_data.append([
                '%s/%s' % (annotation_data.step_id, annotation_data.run_id),
                int(annotation_data.sum['cpu_percent']),
                annotation_data.requested_cores,
                int((annotation_data.sum['rss'] / 1024) / 1024),
                int(annotation_data.total_duration_seconds),
                annotation_data.total_start_time,
                annotation_data.total_end_time
            ])

    if arguments['to_screen']:
        for line in output_data:
            print(line)
    else:
        file_path = arguments['file_path']

        # if file already exists
        if os.path.exists(file_path):
            write_file = False
            while not write_file:
                message = "The given file_path '%s' already exists! " \
                          "Do you want to overwrite this file? [N/y]" % file_path
                response = raw_input(message)

                if response == "N" or response == "":
                    file_path = raw_input("Please give a new file_path:")
                    while file_path == '':
                        file_path = raw_input('Please give a valid file_path:')

                    if os.path.exists(file_path):
                        continue

                    print("New file will be written to '%s'" % file_path)
                    write_file = True
                elif response == 'y':
                    print('File will be overwritten!')
                    write_file = True
                else:
                    print('unknown option')

        write_output_file(output_data, file_path)


if __name__ == '__main__':
    main(args)

class Annotation_Data:

    def __init__(self, file_path):
        file = open(file_path, 'r')
        self.__yaml_data = yaml.load(file, Loader=yaml.FullLoader)
        self.__prepare_data()

    def __prepare_data(self):

        self.main_pid = self.__yaml_data['pid']
        self.run_id = self.__yaml_data['run']['run_id']
        self.step_id = self.__yaml_data['step']['name']

        self.total_start_time = self.__yaml_data['start_time']
        self.total_end_time = self.__yaml_data['end_time']
        self.total_duration = str(self.total_end_time - self.total_start_time)
        self.total_duration_seconds = (self.total_end_time - self.total_start_time).total_seconds()

        self.sum = self.__yaml_data['pipeline_log']['process_watcher']['max']['sum']
        self.requested_cores = self.__yaml_data['step']['cores']

        self.processes = dict()
        for i, process_data in enumerate(self.__yaml_data['pipeline_log']['processes']):
            tmp_data = dict()
            tmp_data['start_time'] = process_data['start_time']
            tmp_data['end_time'] = process_data['end_time']
            tmp_data['duration'] = str(process_data['end_time'] - process_data['start_time'])
            tmp_data['process_pid'] = process_data['pid']

            tmp_data['main_tool'] = 0
            process_name = str(process_data['start_time']) + ' ' + process_data['name']
            if i == 0:
                tmp_data['main_tool'] = 1
                process_name = process_data['name']

            # search in process watcher for information
            if tmp_data['process_pid'] in self.__yaml_data['pipeline_log']['process_watcher']['max']:
                tmp_data.update(self.__yaml_data['pipeline_log']['process_watcher']['max'][tmp_data['process_pid']])

            self.processes[process_name] = dict()
            self.processes[process_name] = tmp_data
