from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class IdentifyAdapters(AbstractStep):
    '''
    Uses AdapterRemoval to identify adapter sequences from paired read data.

    AdapterRemoval (ver. 2.1.7)
    This program searches for and removes remnant adapter sequences from
    your read data.  The program can analyze both single end and paired end
    data.  For detailed explanation of the parameters, please refer to the
    man page.  For comments, suggestions  and feedback please contact Stinus
    Lindgreen (stinus@binf.ku.dk) and Mikkel Schubert (MikkelSch@gmail.com).
    If you use the program, please cite the paper:
    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid
    adapter trimming, identification, and read merging.
    BMC Research Notes, 12;9(1):88.
    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2
    "Pipeline specific "input and output expected to be gzipped"
    '''

    def __init__(self, pipeline):
        super(IdentifyAdapters, self).__init__(pipeline)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

    def runs(self, run_ids_connections_files):
        print('Test')
