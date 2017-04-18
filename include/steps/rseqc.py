import os
import sys
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
        self.add_connection('out/geneBody_coverage.pdf')
        self.add_connection('out/geneBody_coverage_stdout')
        self.add_connection('out/geneBody_coverage_stderr')
        self.add_connection('out/log')

        self.add_connection('out/inner_distance_freq')
        self.add_connection('out/inner_distance_plot')
        self.add_connection('out/inner_distance')
        self.add_connection('out/inner_distance_stdout')
        self.add_connection('out/inner_distance_stderr')

        self.add_connection('out/junction_bed')
        self.add_connection('out/junction_plot')
        self.add_connection('out/junction_xls')
        self.add_connection('out/splice_events')
        self.add_connection('out/splice_junction')
        self.add_connection('out/junction_annotation_stdout')
        self.add_connection('out/junction_annotation_stderr')

        self.add_connection('out/junctionSaturation_pdf')
        self.add_connection('out/junctionSaturation_r')
        self.add_connection('out/junction_saturation_stdout')
        self.add_connection('out/junction_saturation_stderr')

        self.add_connection('out/DupRate_plot_pdf')
        self.add_connection('out/DupRate_plot_r')
        self.add_connection('out/DupRate_pos')
        self.add_connection('out/DupRate_seq')
        self.add_connection('out/DupRate_stdout')
        self.add_connection('out/DupRate_stderr')

        self.add_connection('out/gc_pdf')
        self.add_connection('out/gc_r')
        self.add_connection('out/gc_xls')
        self.add_connection('out/gc_stdout')
        self.add_connection('out/gc_stderr')

        self.require_tool('cat')
        #self.require_tool('cd')
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

    def runs(self, run_ids_connections_files):

        for run_id in run_ids_connections_files.keys():

            alignments = run_ids_connections_files[run_id]['in/alignments']

            with self.declare_run(run_id) as run:
                # todo: can be checked in add_connection('in/') with constraint?
                if len(alignments) != 1:
                    logger.error("Expected exactly one alignment file.")
                    sys.exit(1)

                # todo: why basename and not run_id?
                basename = os.path.basename(alignments[0]).split('.')[0]

                with run.new_exec_group() as exec_group:
                    bam_stat = [
                        self.get_tool('bam_stat.py'),
                        '-i', alignments[0]
                    ]
                    exec_group.add_command(
                        bam_stat,
                        stderr_path=run.add_output_file(
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
                        stdout_path=run.add_output_file(
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
                        stdout_path=run.add_output_file(
                            'read_distribution',
                            basename + '.read_distribution.txt', alignments
                        )
                    )
                out = run.get_output_directory_du_jour_placeholder() + '/' + run_id
                #with run.new_exec_group() as exec_group:
                    # we have to change the working directory because 
                    # geneBody_coverage write a .log file to the current dir
                #    current_dir = os.getcwd()

                #    exec_group.add_command([self.get_tool('cd'), out])
                    #os.chdir(out)

                #    geneBody_coverage = [
                #        self.get_tool('geneBody_coverage.py'),
                #        '-i', alignments[0],
                #        '-r', self.get_option('reference'),
                #        '-o', out
                #    ]
                #    gbc_txt = run_id + '.geneBodyCoverage.txt'
                #    run.add_output_file('geneBody_coverage.txt', gbc_txt, alignments)
                #    gbc_r = run_id + '.geneBodyCoverage.r'
                #    run.add_output_file('geneBody_coverage.r', gbc_r, alignments)
                #    gbc_pdf = run_id + '.geneBodyCoverage.curves.pdf'
                #    run.add_output_file('geneBody_coverage.pdf', gbc_pdf, alignments)
                #    gbc_log = run_id + '.log.txt'
                #    run.add_output_file('log', gbc_log, alignments)

                #    stdout_file = "%s-geneBody_coverage_stdout.txt" % (run_id)
                #    log_stdout = run.add_output_file("geneBody_coverage_stdout",
                #                                     stdout_file, alignments)
                #    stderr_file = "%s-geneBody_coverage_stderr.txt" % (run_id)
                #    log_stderr = run.add_output_file("geneBody_coverage_stderr",
                #                                     stderr_file, alignments)
                    # todo: add log.txt
                #    exec_group.add_command(geneBody_coverage,
                #        stdout_path=log_stdout, stderr_path=log_stderr
                #    )

                    # change dir back to normal working dir
                #    exec_group.add_command([self.get_tool('cd'), current_dir])

                with run.new_exec_group() as exec_group:
                    inner_distance = [
                        self.get_tool('inner_distance.py'),
                        '-i', alignments[0],
                        '-r', self.get_option('reference'),
                        '-o', out
                    ]
                    id_txt = run_id + '.inner_distance.txt'
                    run.add_output_file('inner_distance', id_txt, alignments)
                    id_freq = run_id + '.inner_distance_freq.txt'
                    run.add_output_file('inner_distance_freq', id_freq, alignments)
                    id_plot = run_id + '.inner_distance_plot.r'
                    run.add_output_file('inner_distance_plot', id_plot, alignments)

                    stdout_file = "%s-inner_distance_stdout.txt" % (run_id)
                    log_stdout = run.add_output_file("inner_distance_stdout",
                                                     stdout_file, alignments)
                    stderr_file = "%s-inner_distance_stderr.txt" % (run_id)
                    log_stderr = run.add_output_file("inner_distance_stderr",
                                                     stderr_file, alignments)
                    exec_group.add_command(inner_distance,
                        stdout_path=log_stdout, stderr_path=log_stderr
                    )

                with run.new_exec_group() as exec_group:
                    junction_annotation = [
                        self.get_tool('junction_annotation.py'),
                        '-i', alignments[0],
                        '-r', self.get_option('reference'),
                        '-o', out
                    ]
                    ja_bed = run_id + '.junction.bed'
                    run.add_output_file('junction_bed', ja_bed, alignments)
                    ja_plot = run_id + '.junction_plot.r'
                    run.add_output_file('junction_plot', ja_plot, alignments)
                    ja_xls = run_id + '.junction.xls'
                    run.add_output_file('junction_xls', ja_xls, alignments)
                    ja_events = run_id + '.splice_events.pdf'
                    run.add_output_file('splice_events', ja_events, alignments)
                    ja_junction = run_id + '.splice_junction.pdf'
                    run.add_output_file('splice_junction', ja_junction, alignments)

                    stdout_file = "%s-junction_annotation_stdout.txt" % (run_id)
                    log_stdout = run.add_output_file("junction_annotation_stdout",
                                                     stdout_file, alignments)
                    stderr_file = "%s-junction_annotation_stderr.txt" % (run_id)
                    log_stderr = run.add_output_file("junction_annotation_stderr",
                                                     stderr_file, alignments)

                    exec_group.add_command(junction_annotation,
                        stdout_path=log_stdout, stderr_path=log_stderr
                    )

                with run.new_exec_group() as exec_group:
                    junction_saturation = [
                        self.get_tool('junction_saturation.py'),
                        '-i', alignments[0],
                        '-r', self.get_option('reference'),
                        '-o', out
                    ]
                    js_r = run_id + '.junctionSaturation_plot.r'
                    run.add_output_file('junctionSaturation_r', js_r, alignments)
                    js_pdf = run_id + '.junctionSaturation_plot.pdf'
                    run.add_output_file('junctionSaturation_pdf', js_pdf, alignments)

                    stdout_file = "%s-junction_saturation_stdout.txt" % (run_id)
                    log_stdout = run.add_output_file("junction_saturation_stdout",
                                                     stdout_file, alignments)
                    stderr_file = "%s-junction_saturation_stderr.txt" % (run_id)
                    log_stderr = run.add_output_file("junction_saturation_stderr",
                                                     stderr_file, alignments)

                    exec_group.add_command(junction_saturation,
                        stdout_path=log_stdout, stderr_path=log_stderr
                    )

                with run.new_exec_group() as exec_group:
                    read_duplication = [
                        self.get_tool('read_duplication.py'),
                        '-i', alignments[0],
                        '-o', out
                    ]
                    rd_pdf = run_id + '.DupRate_plot.pdf'
                    run.add_output_file('DupRate_plot_pdf', rd_pdf, alignments)
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
                        stdout_path=log_stdout, stderr_path=log_stderr
                    )

                with run.new_exec_group() as exec_group:
                    read_gc = [
                        self.get_tool('read_GC.py'),
                        '-i', alignments[0],
                        '-o', out
                    ]
                    gc_pdf = run_id + '.GC_plot.pdf'
                    run.add_output_file('gc_pdf', gc_pdf, alignments)
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
                        stdout_path=log_stdout, stderr_path=log_stderr
                    )
