from logging import getLogger
from uaperrors import UAPError

logger = getLogger("uap_logger")

class ConnectionsCollector(object):
    def __init__(self, step_name=None):
        self.step_name = step_name
        self.connections = dict()
        self._current_run_id = None
        self.used_current_run_id = False
        self._by_cons_empty = dict()
        self._by_cons_none_empty = dict()
        self._con_of_all_runs = None
        self.existing_connections = set()

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
            raise UAPError('No run id given.')
        return run_id

    def add_empty(self, connection, run_id=None):
        run_id = self.init_run_id(run_id)
        self.connections[run_id][connection] = [None]
        self._by_cons_empty.setdefault(connection, set())
        self._by_cons_empty[connection].add(un_id)
        self.existing_connections.add(connection)
        self._con_of_all_runs = None # reset cache
        logger.debug("Found connection %s which is declared empty" %
                     (connection))

    def add_connection(self, connection, files, run_id=None):
        run_id = self.init_run_id(run_id)
        self.connections[run_id].setdefault(connection, list())
        self.connections[run_id][connection].extend(files)
        self._by_cons_none_empty.setdefault(connection, set())
        self._by_cons_none_empty[connection].add(run_id)
        self.existing_connections.add(connection)
        self._con_of_all_runs = None # reset cache
        logger.debug("Found %s to connect to %s in run %s." %
                (self.step_name, connection, run_id))

    def get_connection(self, connection, run_id=None):
        run_id = self.init_run_id(run_id)
        cons = self.connections[run_id]
        if connection not in cons.keys():
            raise UAPError('The input run %s of %s has no connection %s.' %
                    (run_id, self.step_name, connection))
        return cons[run_id]

    def _runs_with_connection(self, connection, with_empty=True):
        runs = set()
        if connection in self._by_cons_none_empty.keys():
            runs = runs.union(self._by_cons_none_empty[connection])
        if with_empty is True:
            if connection in self._by_cons_empty.keys():
                runs = runs.union(self._by_cons_empty[connection])
        return runs

    def get_runs_with_connections(self, connections, with_empty=True):
        '''
        Returns all run ids that have requested all connections.
        '''
        if isinstance(connections, str):
            return self._runs_with_connection(connections, with_empty)
        cons = list(connections)
        run_ids = self._runs_with_connection(connections[0], with_empty)
        for con in connections[1:]:
            run_ids.intersection(self._runs_with_connection(con, with_empty))
        return run_ids

    def look_for_unique(self, connection, default=None):
        '''
        Looks for a unique file in the connection and returns it.
        E.g., to find a reference assembly among all parent runs.
        If all runs come with the file and no default is set it returns None.
        If more then one but not all runs come with the connection an UAPError is raised.
        '''
        if self.connection_exists(connection) and default is not None:
            raise UAPError('In step %s runs come with %s but it is set '
                'to %s in the config.' % (self.step_name, connection, default))

        if self.all_runs_have_connection(connection):
            return default

        ref_run = self.get_runs_with_connections(connection, with_empty=False)
        if len(ref_run) > 1:
            UAPError('More then one but not all runs come with %s.' % connection)
        elif len(ref_run) == 1:
            if default is not None:
                raise UAPError('In step %s, value supplied by connection %s but'
                        'option is set to %s.' % (self.step_name, connection, default))
            ref_run = ref_run.pop()
            con_value = self.connections[ref_run][connection]
            if len(con_value) > 1:
                raise UAPError('In step %s more than one file is passed with %s.' %
                        (self.step_name, connection))
            return con_value[0]
        return default

    def connection_exists(self, connection):
        return connection in self.existing_connections

    def all_runs_have_connection(self, connection):
        '''
        Returns True or False.
        '''
        if self._con_of_all_runs is None:
            # calculate and cache connections of all runs
            con_list = self.connections.values()
            if len(con_list) > 1:
                self._con_of_all_runs = set(con_list[0]).intersection(*con_list)
            elif len(con_list) == 1:
                self._con_of_all_runs = set(con_list[0])
            else:
                self._con_of_all_runs = set()
        return connection in self._con_of_all_runs

    def add_default_ins(self, out_connection, files):
        in_connection = out_connection.replace('out/', 'in/')
        self.add_connection(in_connection, files)

    def __getitem__(self, run_id):
        if run_id not in self.connections.keys():
            raise KeyError('In step %s there is no connection for run %s.' %
                    (self.step_name, run_id))
        return self.connections[run_id]

    def keys(self):
        return self.connections.keys()

    def values(self):
        return self.connections.values()

    def items(self):
        return self.connections.items()
