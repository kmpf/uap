from logging import getLogger
from uaperrors import UAPError

logger = getLogger("uap_logger")


class ConnectionsCollector(object):
    '''
    A ConnectionsCollector helps to collect and query file connections.
    An instance `cc` is generated for each step and passed to its
    ``runs`` method. For backwards compatibility reasons it can be
    queryed like a dictionary ``cc[run_id][connection]``. ``cc``
    can be used in the course of a step to dicide how to use the
    input runs and connections.
    '''

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
        '''
        Set a default run id to use in subsequent calls.
        '''
        self._init_run_id(run_id)
        self._current_run_id = run_id
        self.used_current_run_id = False

    def _init_run_id(self, run_id=None):
        if run_id is None:
            run_id = self._current_run_id
        else:
            self.connections.setdefault(run_id, dict())
        if run_id is None:
            raise UAPError('No run id given.')
        return run_id

    def add_empty(self, connection, run_id=None):
        '''
        Makes a connection withou files.
        '''
        run_id = self._init_run_id(run_id)
        self.connections[run_id][connection] = [None]
        self._by_cons_empty.setdefault(connection, set())
        self._by_cons_empty[connection].add(run_id)
        self.existing_connections.add(connection)
        self._con_of_all_runs = None  # reset cache
        logger.debug("Found connection %s which is declared empty" %
                     (connection))

    def add_connection(self, connection, files, run_id=None):
        '''
        Saves the names in ``files`` for a new ``connection``.
        '''
        if not isinstance(files, list):
            raise UAPError('"files" must be a list but is a %s' %
                           files.__class__.__name__)
        run_id = self._init_run_id(run_id)
        if not isinstance(connection, str):
            raise UAPError('The passed connection must be a string.')
        if not connection.startswith('in/'):
            raise UAPError('Input connections muss start with "in/".')
        self.connections[run_id].setdefault(connection, list())
        self.connections[run_id][connection].extend(files)
        self._by_cons_none_empty.setdefault(connection, set())
        self._by_cons_none_empty[connection].add(run_id)
        self.existing_connections.add(connection)
        self._con_of_all_runs = None  # reset cache
        logger.debug("Found %s to connect %s with run %s." %
                     (self.step_name, connection, run_id))

    def connect(self, parent, child, connections=None):
        '''
        Makes connections between parent and child step and returns
        the utilized parent connections as "parent/name".
        The passed connections need to be in the format as in the
        pipline configuration file. If no connections are passed it
        connects all equally named connections.
        '''
        parent_name = parent.get_step_name()
        parent_out_conns = parent.get_out_connections(strip_prefix=True)
        child_name = child.get_step_name()
        must_connect = set()
        if connections is None:
            # get equally named connections
            ins = child.get_in_connections(strip_prefix=True)
            conns = parent_out_conns.intersection(ins)
            if len(conns) == 0:
                logger.warning(
                    'There are no default connections between '
                    '%s and its dependency %s. The parent out connections '
                    'are %s and the child in connections are %s.' %
                    (parent_name, child_name, list(parent_out_conns), list(ins)))
                return 0
            make_connections = [
                ('in/%s' % conn,
                 'out/%s' % conn,
                 '%s/%s' % (parent_name, conn)) for conn in conns]
        elif isinstance(connections, dict):
            # extract connections from config
            make_connections = list()
            pre_len = len(parent_name) + 1
            for in_conn, out_conns in connections.items():
                if not isinstance(out_conns, list):
                    out_conns = [out_conns]
                for out_conn in out_conns:
                    if out_conn.startswith(parent_name + '/'):
                        stripped_out_conn = out_conn[pre_len:]
                        p_out = 'out/%s' % stripped_out_conn
                        if stripped_out_conn not in parent_out_conns:
                            avail = list(parent_out_conns)
                            raise UAPError(
                                'The connection "%s" set in "%s" '
                                'is not an out connection of "%s". '
                                'Available out connections are: %s' %
                                (out_conn, child_name, parent_name, avail))
                        make_connections.append((in_conn, p_out, out_conn))
                        must_connect.add(out_conn)
        else:
            raise UAPError('The passed connections need to be a dictionay.')

        # make the connections
        used_conns = set()
        for parent_run_id in parent.get_runs():
            self.switch_run_id(parent_run_id)
            parent_run = parent.get_run(parent_run_id)
            for in_conn, out_conn, parent_con in make_connections:
                if out_conn not in parent_run.get_out_connections():
                    continue
                output_files = parent_run\
                    .get_output_files_abspath_for_out_connection(out_conn)
                self.add_connection(in_conn, output_files)
                used_conns.add(parent_con)

        missing = must_connect - used_conns
        if missing:
            raise UAPError(
                'The connection(s) %s required by the step "%s" are '
                'optional output of step "%s" and not produced.' %
                (list(missing), child_name, parent_name))

        return used_conns

    def get_connection(self, connection, run_id=None):
        '''
        Returns a list of file names for the ``connection``.
        '''
        run_id = self._init_run_id(run_id)
        cons = self.connections[run_id]
        if connection not in cons.keys():
            raise UAPError('The input run %s of %s has no connection %s.' %
                           (run_id, self.step_name, connection))
        return cons[run_id]

    def connection_items(self, connection):
        '''
        Returns all (run_id, [files]) pairs for the given connection.
        '''
        for run_id, conns in self.items():
            if connection in conns:
                yield run_id, conns[connection]

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
        run_ids = self._runs_with_connection(cons[0], with_empty)
        for con in cons[1:]:
            runs_of_con = self._runs_with_connection(con, with_empty)
            run_ids = run_ids.intersection(runs_of_con)
        return run_ids

    def get_runs_with_any(self, connections, with_empty=True):
        '''
        Returns all run ids with any of the given connections.
        '''
        if isinstance(connections, str):
            connections = [connections]
        cons = list(connections)
        run_ids = self._runs_with_connection(cons[0], with_empty)
        for con in cons[1:]:
            runs_of_con = self._runs_with_connection(con, with_empty)
            run_ids = run_ids.union(runs_of_con)
        return run_ids

    def get_runs_without_any(self, connections, with_empty=True):
        '''
        Returns runs that have none of the passed connections.
        '''
        runs = set(self.connections.keys())
        runs_with_any = self.get_runs_with_any(connections)
        return runs.difference(runs_with_any)

    def look_for_unique(self, connection, default=None):
        '''
        Looks for a unique file in the connection and returns it.
        E.g., to find a reference assembly among all parent runs.
        If all runs come with the file and no default is set it returns None.
        If more then one but not all runs come with the connection an UAPError is raised.
        '''
        if self.connection_exists(connection) and default is not None:
            raise UAPError(
                'In step %s runs come with %s but it is set '
                'to %s in the config.' %
                (self.step_name, connection, default))

        if self.all_runs_have_connection(connection):
            return default

        ref_run = self.get_runs_with_connections(connection, with_empty=False)
        if len(ref_run) > 1:
            UAPError(
                'More then one but not all runs come with %s.' %
                connection)
        elif len(ref_run) == 1:
            if default is not None:
                raise UAPError(
                    'In step %s, value supplied by connection %s but'
                    'option is set to %s.' %
                    (self.step_name, connection, default))
            ref_run = ref_run.pop()
            con_value = self.connections[ref_run][connection]
            if len(con_value) > 1:
                raise UAPError(
                    'In step %s more than one file is passed with %s.' %
                    (self.step_name, connection))
            return con_value[0]
        return default

    def connection_exists(self, connection):
        '''
        Returns a logical indication whether the requested connection exists.
        '''
        return connection in self.existing_connections

    def all_runs_have_connection(self, connection):
        '''
        Returns a logical indication whether all saved runs have the queried
        ``connection``.
        '''
        if self._con_of_all_runs is None:
            # calculate and cache connections of all runs
            con_list = self.connections.values()
            self._con_of_all_runs = set.intersection(
                *(set(val) for val in con_list))
        if isinstance(connection, list):
            return all(con in self._con_of_all_runs for con in connection)
        return connection in self._con_of_all_runs

    def __getitem__(self, run_id):
        if run_id not in self.connections.keys():
            raise KeyError('In step %s there is no connection for run %s.' %
                           (self.step_name, run_id))
        return self.connections[run_id]

    def keys(self):
        return self.connections.keys()

    def values(self):
        '''
        Emulates dict.values().
        '''
        return self.connections.values()

    def items(self):
        '''
        Emulates dict.items().
        '''
        return iter(sorted(self.connections.items()))
