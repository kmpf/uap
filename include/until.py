from logging import getLogger
from uaperrors import UAPError

logger = getLogger("uap_logger")

class ConnectionsCollector(Object):
    def __init__(self, name=None):
        self.name = name
        self.connections = dict()
        self._current_run_id = None
        self.used_current_run_id = False

    def switch_run_id(self, run_id):
        self.init_run_id(run_id)
        self._current_run_id = run_id
        self.used_current_run_id = False

    def init_run_id(self, run_id=None):
        if run_id is None:
            run_id = self._current_run_id
        else:
            self.connections.setdefault(run_id, dict())
        if run_id is None:
            UAPError('No run id given.')
        return run_id

    def add_empty(self, connection, run_id=None):
        run_id = self.init_run_id(run_id)
        self.connections[run_id][connection] = [None]
        logger.debug("Found connection %s which is declared empty" %
                     (connection))

    def add_connection(self, connection, files, run_id=None):
        run_id = self.init_run_id(run_id)
        logger.debug("Found %s to connect to %s in run %s." %
                (self.name, connection, run_id))
        self.connections[run_id].setdefault(connection, list())
        self.connections[run_id].[connection].extend(files)

    def add_default_ins(self, out_connection, files):
        in_connection = out_connection.replace('out/', 'in/')
        self.add_connection(in_connection, files)
