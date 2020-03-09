import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class RSeQC(AbstractStep):
    '''
    The RSeQC step can be used to evaluate aligned reads in a BAM file. RSeQC
    does not only report raw sequence-based metrics, but also quality control
    metrics like read distribution, gene coverage, and sequencing depth.
    '''

    def __init__(self, pipeline):
        super(RSeQC, self).__init__(pipeline)

        self.set_cores(4)

        self.add_connection('in/alignments')

        self.add_connection('out/bam_stat')
        self.add_connection('out/infer_experiment')
        self.add_connection('out/read_distribution')

        self.add_connection('out/geneBody_coverage.txt')
        self.add_connection('out/geneBody_coverage.r')
        self.add_connection('out/geneBody_coverage_stdout')
        self.add_connection('out/geneBody_coverage_stderr')

        self.add_connection('out/inner_distance_freq', optional=True)
        self.add_connection('out/inner_distance_plot', optional=True)
        self.add_connection('out/inner_distance', optional=True)
        self.add_connection('out/inner_distance_stdout', optional=True)
        self.add_connection('out/inner_distance_stderr', optional=True)

        self.add_connection('out/junction_bed')
        self.add_connection('out/junction_plot')
        self.add_connection('out/junction_xls')

        self.add_connection('out/junction_annotation_stdout')
        self.add_connection('out/junction_annotation_stderr')

        self.add_connection('out/junctionSaturation_r')
        self.add_connection('out/junction_saturation_stdout')
        self.add_connection('out/junction_saturation_stderr')

        self.add_connection('out/DupRate_plot_r')
        self.add_connection('out/DupRate_pos')
        self.add_connection('out/DupRate_seq')
        self.add_connection('out/DupRate_stdout')
        self.add_connection('out/DupRate_stderr')

        self.add_connection('out/gc_r')
        self.add_connection('out/gc_xls')
        self.add_connection('out/gc_stdout')
        self.add_connection('out/gc_stderr')

        self.require_tool('cat')
        self.require_tool('rm')
        self.require_tool('bam_stat.py')
        self.require_tool('infer_experiment.py')
        self.require_tool('read_distribution.py')
        self.require_tool('geneBody_coverage.py')
        self.require_tool('inner_distance.py')
        self.require_tool('junction_annotation.py')
        self.require_tool('junction_saturation.py')
        self.require_tool('read_duplication.py')
        self.require_tool('read_GC.py')

        self.add_option('reference', str, optional=False,
                        description="Reference gene model in bed fomat. "
                        "[required]")

        self.add_option('treatAs', str, optional=False,
                        choices = ["single", "paired"],
                        description="Some modules in rseqc  need paired end data"
                        "an fail otherwise on single end [required]")

    def runs(self, run_ids_connections_files):

        for run_id in run_ids_connections_files.keys():

            alignments = run_ids_connections_files[run_id]['in/alignments']

            with self.declare_run(run_id) as run:
                with run.new_exec_group() as exec_group:
                    bam_stat = [
                        self.get_tool('bam_stat.py'),
                        '-i', alignments[0]
                    ]
                    exec_group.add_command(
                        bam_stat,
                        stderr_path=run.add_output_file(
                            'bam_stat', run_id + '.bam_stats.txt', alignments
                        )
                    )

                with run.new_exec_group() as exec_group:
                    infer_experiment = [
                        self.get_tool('infer_experiment.py'),
                        '-i', alignments[0],
                        '-r', os.path.abspath(self.get_option('reference'))
                    ]
                    exec_group.add_command(
                        infer_experiment,
                        stdout_path=run.add_output_file(
                            'infer_experiment',
                            run_id + '.infer_experiment.txt', alignments
                        )
                    )

                with run.new_exec_group() as exec_group:
                    read_distribution = [
                        self.get_tool('read_distribution.py'),
                        '-i', alignments[0],
                        '-r', os.path.abspath(self.get_option('reference'))
                    ]
                    exec_group.add_command(
                        read_distribution,
                        stdout_path=run.add_output_file(
                            'read_distribution',
                            run_id + '.read_distribution.txt', alignments
                        )
                    )

                with run.new_exec_group() as exec_group:
                    current_dir = os.getcwd()

                    geneBody_coverage = [
                        self.get_tool('geneBody_coverage.py'),
                        '-i', alignments[0],
                        '-r', os.path.abspath(self.get_option('reference')),
                        '-o', run_id
                    ]
                    gbc_txt = run_id + '.geneBodyCoverage.txt'
                    run.add_output_file('geneBody_coverage.txt',
                                        gbc_txt, alignments)
                    gbc_r = run_id + '.geneBodyCoverage.r'
                    run.add_output_file('geneBody_coverage.r',
                                        gbc_r, alignments)

                    stdout_file = "%s-geneBody_coverage_stdout.txt" % (run_id)
                    log_stdout = run.add_output_file(
                        "geneBody_coverage_stdout", stdout_file, alignments)
                    stderr_file = "%s-geneBody_coverage_stderr.txt" % (run_id)
                    log_stderr = run.add_output_file(
                        "geneBody_coverage_stderr", stderr_file, alignments)

                    exec_group.add_command(geneBody_coverage,
                                           stdout_path=log_stdout,
                                           stderr_path=log_stderr)


                    if self.get_option('treatAs') == 'paired':
                        inner_distance = [
                            self.get_tool('inner_distance.py'),
                            '-i', alignments[0],
                            '-r', os.path.abspath(self.get_option('reference')),
                            '-o', run_id
                        ]
                        id_txt = run_id + '.inner_distance.txt'
                        run.add_output_file('inner_distance', id_txt, alignments)
                        id_freq = run_id + '.inner_distance_freq.txt'
                        run.add_output_file('inner_distance_freq',
                                            id_freq, alignments)
                        id_plot = run_id + '.inner_distance_plot.r'
                        run.add_output_file('inner_distance_plot',
                                            id_plot, alignments)

                        stdout_file = "%s-inner_distance_stdout.txt" % (run_id)
                        log_stdout = run.add_output_file("inner_distance_stdout",
                                                     stdout_file, alignments)

                        stderr_file = "%s-inner_distance_stderr.txt" % (run_id)
                        log_stderr = run.add_output_file("inner_distance_stderr",
                                                     stderr_file, alignments)

                        exec_group.add_command(inner_distance,
                                               stdout_path=log_stdout,
                                               stderr_path=log_stderr)

                with run.new_exec_group() as exec_group:
                    junction_annotation = [
                        self.get_tool('junction_annotation.py'),
                        '-i', alignments[0],
                        '-r', os.path.abspath(self.get_option('reference')),
                        '-o', run_id
                    ]
                    ja_bed = run_id + '.junction.bed'
                    run.add_output_file('junction_bed', ja_bed, alignments)
                    ja_plot = run_id + '.junction_plot.r'
                    run.add_output_file('junction_plot', ja_plot, alignments)
                    ja_xls = run_id + '.junction.xls'
                    run.add_output_file('junction_xls', ja_xls, alignments)

                    stdout_file = "%s-junction_annotation_stdout.txt" \
                        % (run_id)
                    log_stdout = run.add_output_file(
                        "junction_annotation_stdout", stdout_file, alignments)
                    stderr_file = "%s-junction_annotation_stderr.txt" \
                        % (run_id)
                    log_stderr = run.add_output_file(
                        "junction_annotation_stderr", stderr_file, alignments)

                    exec_group.add_command(junction_annotation,
                                           stdout_path=log_stdout,
                                           stderr_path=log_stderr)

                with run.new_exec_group() as exec_group:
                    junction_saturation = [
                        self.get_tool('junction_saturation.py'),
                        '-i', alignments[0],
                        '-r', os.path.abspath(self.get_option('reference')),
                        '-o', run_id
                    ]
                    js_r = run_id + '.junctionSaturation_plot.r'
                    run.add_output_file('junctionSaturation_r',
                                        js_r, alignments)

                    stdout_file = "%s-junction_saturation_stdout.txt" \
                        % (run_id)
                    log_stdout = run.add_output_file(
                        "junction_saturation_stdout", stdout_file, alignments)
                    stderr_file = "%s-junction_saturation_stderr.txt" \
                        % (run_id)
                    log_stderr = run.add_output_file(
                        "junction_saturation_stderr", stderr_file, alignments)

                    exec_group.add_command(junction_saturation,
                                           stdout_path=log_stdout,
                                           stderr_path=log_stderr)

                with run.new_exec_group() as exec_group:
                    read_duplication = [
                        self.get_tool('read_duplication.py'),
                        '-i', alignments[0],
                        '-o', run_id
                    ]
                    rd_r = run_id + '.DupRate_plot.r'
                    run.add_output_file('DupRate_plot_r', rd_r, alignments)
                    rd_pos = run_id + '.pos.DupRate.xls'
                    run.add_output_file('DupRate_pos', rd_pos, alignments)
                    rd_seq = run_id + '.seq.DupRate.xls'
                    run.add_output_file('DupRate_seq', rd_seq, alignments)

                    stdout_file = "%s-read_duplication_stdout.txt" % (run_id)
                    log_stdout = run.add_output_file("DupRate_stdout",
                                                     stdout_file, alignments)
                    stderr_file = "%s-read_duplication_stderr.txt" % (run_id)
                    log_stderr = run.add_output_file("DupRate_stderr",
                                                     stderr_file, alignments)

                    exec_group.add_command(read_duplication,
                                           stdout_path=log_stdout,
                                           stderr_path=log_stderr)

                with run.new_exec_group() as exec_group:
                    read_gc = [
                        self.get_tool('read_GC.py'),
                        '-i', alignments[0],
                        '-o', run_id
                    ]
                    gc_r = run_id + '.GC_plot.r'
                    run.add_output_file('gc_r', gc_r, alignments)
                    gc_xls = run_id + '.GC.xls'
                    run.add_output_file('gc_xls', gc_xls, alignments)

                    stdout_file = "%s-read_gc_stdout.txt" % (run_id)
                    log_stdout = run.add_output_file("gc_stdout",
                                                     stdout_file, alignments)
                    stderr_file = "%s-read_gc_stderr.txt" % (run_id)
                    log_stderr = run.add_output_file("gc_stderr",
                                                     stderr_file, alignments)

                    exec_group.add_command(read_gc,
                                           stdout_path=log_stdout,
                                           stderr_path=log_stderr)

