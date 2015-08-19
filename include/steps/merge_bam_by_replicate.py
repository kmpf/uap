import sys
import yaml

from ..abstract_step import *
from .. import process_pool

class Merge_Bam_By_Replicate (AbstractStep):

    def __init__(self, pipeline):
        super(Merge_Bam_By_Replicate, self).__init__(pipeline)

        self.set_cores (1)

        self.add_connection('in/alignments')
        self.add_connection('out/alignments')
        self.add_connection('out/indices')
        #self.add_connection('out/log')

        self.add_option('replicate_groups', dict)

        self.require_tool('bamtools')
        self.require_tool('samtools')
        self.require_tool('cat')
        ## ? also sort  pigz and the like ??


    def declare_runs (self):


        ## Check whether replicate grouping is unique, i.e. whether the given "member names" match only one
        ## of the input alignments
        replicate_groups=self.get_option('replicate_groups')
        ### Comment when productive
        ##pp = pprint.PrettyPrinter(indent=4)
        #print reduce(lambda x, y: x+y,replicate_groups.values())
        ##pp.pprint(replicate_groups)
        ###pp.pprint(list(self.get_run_ids_and_input_files_for_connection ('in/alignments')))
        ##pp.pprint (replicate_groups.values())
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection ('in/alignments'):
            #print run_id
            ## check whether run_id is matched only once by replicate group members
            match_list=map(lambda a: a in run_id, reduce(lambda x, y: x+y,replicate_groups.values())) ##reduce used to  flatten the list of lists
            if (match_list.count(True) > 1):
                pp.pprint(match_list)
                raise StandardError("Ambiguous replicate group members: run_id matched more than once %s: %s" % (run_id, self))

            

            ## repeat for paths
            match_list=map(lambda a: a in input_paths, reduce(lambda x, y: x+y,replicate_groups.values()))
            if (match_list.count(True) > 1):
                raise StandardError("Ambiguous replicate group members: input_paths matched more than once %s: %s" % (run_id, self))


        ## Setup runs only for the replicate groups instead of all input files
        ##run_ids, input_paths=self.get_run_ids_and_input_files_for_connection ('in/alignments')
            
        result=self.get_input_run_info_for_connection ('in/alignments')
        #print "Runkeys:\n"
        #print result['runs'].keys()
        #print result['runs']['T-Cell_BC-16_24h-R1'].values()
        #print result['runs'].items()

        run_ids=result['runs'].keys()
        #print run_ids
        
        for replicate_group in replicate_groups.keys():
            with self.declare_run (replicate_group) as run:
                run.new_exec_group()
                input_run_ids_for_replicate_group = list()
                input_file_paths_for_replicate_g = list()
                ## Find run_ids that match the replicate group members and the corresponding input files
                for r in run_ids:
                    for g in replicate_groups[replicate_group]:
                        if (g in r):
                            input_run_ids_for_replicate_group.append(r)
                            input_file_paths_for_replicate_g.append (reduce(lambda x, y: x+y, result['runs'][r].values() ))
                            
                            
                            ## Flatten list again - not clear to me why necessary repeatedly
                            input_file_paths_for_replicate_group = reduce(lambda x, y: x+y, input_file_paths_for_replicate_g)
                            
                            
                            #print replicate_group
                            #print input_run_ids_for_replicate_group
                            #print input_file_paths_for_replicate_group
                            #print input_file_paths_for_replicate_group
            
                run.add_output_file ('alignments', "%s-merged.bam" % replicate_group, input_file_paths_for_replicate_group)
                run.add_output_file ('indices', "%s-merged.bam.bai" % replicate_group, input_file_paths_for_replicate_group)
                #run.add_output_file ('log', "%s-merged-log.txt" % replicate_group, input_file_paths_for_replicate_group)
                
            
    def execute (self, run_id, run):
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                ## Obtain the output files
                sorted_bam_path = run.get_single_output_file_for_annotation('alignments')
                sorted_bai_path = run.get_single_output_file_for_annotation('indices')
                unsorted_bam_path = self.get_temporary_path('merge_bam_unsorted', 'output')

                ## Merge bam files using bamtools merge - maybe replace by more performant solutin using cat
                bamtools =[self.get_tool('bamtools'), 'merge']
                bamtools.extend (['-forceCompression'])
                ## Append the input files, which are all input file linked to the output file
                for file in run.get_input_files_for_output_file(sorted_bam_path):
                    bamtools.extend (['-in', file])
                    
                #bamtools.extend ( map (lambda a: '-in ' +a, run.get_input_files_for_output_file(sorted_bam_path)))
                #bamtools.extend ([ '-out ', unsorted_bam_path])
                
                pipeline.append (bamtools, stdout_path=unsorted_bam_path)#, hints={'writes': [unsorted_bam_path]})
                
        ## Sort the bam file by chromosome
        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:
                cat = [self.get_tool('cat'), unsorted_bam_path]
                samtools = [self.get_tool('samtools'), 'sort']
                samtools.extend(['-', sorted_bam_path[:-4]])
                
                pipeline.append(cat)
                pipeline.append(samtools, hints = {'writes': [sorted_bam_path]})

        # samtools index
        
        with process_pool.ProcessPool(self) as pool:
            pool.launch([self.get_tool('samtools'), 'index', sorted_bam_path, '/dev/stdout'],
                stdout_path = sorted_bai_path, hints = {'reads': [sorted_bam_path]})

        if os.path.exists(unsorted_bam_path):
            os.unlink(unsorted_bam_path)
            
        
