import sys
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class RSeQC(AbstractStep):
    '''
    The RSeQC step can be used to evaluate aligned reads in a BAM file. RSeQC
    does not only report raw sequence-based metrics, but also quality control
    metrics like read distribution, gene coverage, and sequencing depth.
    '''

    def __init__(self, pipeline):
        super(RSeQC, self).__init__(pipeline)
        
        self.set_cores(1)
        
        self.add_connection('in/alignments')
        self.add_connection('out/bam_stat')
	self.add_connection('out/infer_experiment')
	self.add_connection('out/read_distribution')
                
        self.require_tool('cat')
        self.require_tool('bam_stat.py')
	self.require_tool('infer_experiment.py')
	self.require_tool('read_distribution.py')

	self.add_option('reference', str, optional=False, description=
                        "Reference gene model in bed fomat. [required]"
        )
        
    def runs(self, run_ids_connections_files):
        
        for run_id  in run_ids_connections_files.keys():

            alignments = run_ids_connections_files[run_id]['in/alignments']

            with self.declare_run(run_id) as run:

                if len(alignments) != 1:
                    logger.error("Expected exactly one alignment file.")
                    sys.exit(1)

                basename = os.path.basename(alignments[0]).split('.')[0]

                with run.new_exec_group() as exec_group:
                    bam_stat = [
                        self.get_tool('bam_stat.py'),
                        '-i', alignments[0]
                    ]
                    exec_group.add_command(
                        bam_stat,
                        stderr_path = run.add_output_file(
                            'bam_stat', basename + '.bam_stats.txt', alignments
                        )
                    )

                with run.new_exec_group() as exec_group:
                    infer_experiment = [
                        self.get_tool('infer_experiment.py'),
                        '-i', alignments[0],
                        '-r', self.get_option('reference')
                    ]
                    exec_group.add_command(
                        infer_experiment,
                        stdout_path = run.add_output_file(
                            'infer_experiment',
                            basename + '.infer_experiment.txt', alignments
                        )
                    )

                with run.new_exec_group() as exec_group:
                    read_distribution = [
                        self.get_tool('read_distribution.py'),
                        '-i', alignments[0],
                        '-r', self.get_option('reference')
                    ]
                    exec_group.add_command(
                        read_distribution,
                        stdout_path = run.add_output_file(
                            'read_distribution',
                            basename + '.read_distribution.txt', alignments
                        )
                    )
