__author__ = 'Jonas I Liechti'
import numpy.random as random
from numpy import array
import numpy as np
from Queue import Queue
from copy import copy

# Set of possible distributions for the degree.
Distribution = {
    'poisson': random.poisson,
    'normal': random.normal,
    'binomial': random.binomial,
    'exponential': random.exponential,
    'geometric': random.geometric,
    'gamma': random.gamma,
    'power': random.power,
    'weibull': random.weibull,
    'negative_binomial': random.negative_binomial
}

allowed_dists = Distribution.keys()


class InvalidArgumentError(Exception):
    def __init__(self, msg):
        self.msg = msg


class NotImplementedError(Exception):
    def __init__(self, msg):
        self.msg = msg


class Node():
    all_ = list()

    class NoNodeError(Exception):
        pass

    def __init__(
            self, uid=None, contacts=None, neighbours=None, start=None,
            stop=None, prop=None, get_from_prop=None
            ):
        """
        This class defines a single node.

        Arguments:
            - uid: unique identifier of the node
            - contacts: a list of contacts the node has.
            - neighbours: a list of neighbours of the node
        :param uid: unique id for each host in a contact_structure
        :param contacts: A list of either lists/tuples or dicts. Each element
            describes a contact that must have a partner node (use an_id
            attribute) and can have a duration and/or start stop and further
            info.
        :param neighbours: A list of either host uid's indicating all
            the neighbours of the node.
        :param prop: Optional argument. If provided, must be a dict.
        :param get_from_prop: List of strings specifying which keys from the
            dict provided in prop should be made attributes of the node.
        :return:
        """
        if self not in self.__class__.all_:
            self.__class__.all_.append(self)
            self._id = len(self.__class__.all_) - 1
        else:
            self_id = self.__class__.all_.index(self)
        if uid is not None:
            self.uid = uid
        else:
            self.uid = _id
        self.neighbours = neighbours
        self.contacts = []
        if contacts is not None:
            self.add_contacts(contacts)
        self.start=start
        self.stop=stop
        if prop is not None:
            if get_from_prop is not None:
                for g_prop in get_from_prop:
                    setattr(self, g_prop, prop.pop(g_prop, None))

    def __getstate__(self):
        instance_dict = self.__dict__
        return instance_dict

    def __setstate__(self, input_dict):
        self.__dict__ = input_dict
        if self not in self.__class__.all_:
            self.__class__.all_.append(self)
            self._id = len(self.__class__.all_) - 1
        else:
            self_id = self.__class__.all_.index(self)

    @classmethod
    def _no_instances(cls):
        if not len(cls.all_):
            raise cls.NoNodeError(
                    'The set of Nodes is empty. Add a node first'
                    )
        else:
            return 1

    @classmethod
    def get_node(cls, uid):
        try:
            cls._no_instances()
        except cls.NoNodeError:
            return None
        for node in cls.all_:
            if uid == node.uid:
                return node
        return None

    @staticmethod
    def add_events(to_add, event):
        for a_key, val in event.iteritems():
            try:
                to_add[a_key].append(val)
            except KeyError:
                to_add[a_key] = [val]

    def add_contacts(self, contacts):
        self.add_events(self.contacts, contacts)

    #def add_manipulations(self, manipulations):
    #    self.add_events(self.manipulations, manipulations)


class _Graph():
    """
        This is the basic class for a graph, containing but
        :param nodes:
        :param edges: Is either a set of tuples containing node ids or node ids
             and start and stop times
        :return:
    """
    def __init__(self, nodes=None, n=None, edges=None, degrees=None):
        self._nodes = nodes if nodes is not None else np.array([])
        self._edges = edges if edges is not None else np.array([])
        self.degrees = None
        if degrees is not None:
            self.degrees = degrees
        # to do: this is not ideal when working with np.array's
        if edges is not None:
            nn = [[] for _ in xrange(len(self._nodes))]
            self.degrees = []
            for edge in self._edges:
                nn[edge[0]].append(edge[1])
                nn[edge[1]].append(edge[0])
            for i in xrange(len(nn)):
                self.degrees.append(len(list(set(nn[i]))))
            del nn
        self.n = None
        if n is not None:
            self.n = n
        if np.size(self._nodes):
            self.n = len(self._nodes)
        self.info = {}  # dict containing some info


class TemporalGraph(_Graph):
    """
        Source is the actual data we want to import. In params we can specify
            if we just want to look at part of the data (e.g. provide 't_start'
            and 't_stop').
        :param nodes: Is a list of instances of the Node class.
            Note: you should pass Node.nodes here. If this argument is not
            provided either the events or the source arguments need to be
            provided.
        :param node_import: A list of node attributes that will be converted to
            an attribute of self (in form of a list).
        :param events: a list of events, each element must contain a start,
            stop node1 and node2. The elements can be lists itself or dicts.
            If events is not provided the host or temporal_graph attribute must
            be provided.
        :param event_keys: A dictionary mapping the following keys:
            'start', 'stop', 'node1', 'node2'. This attribute must be provided
            if the events attribute is not None. The corresponding values to
            these keys must allow to extract the content from each individual
            event from the events attribute. So if events is a list of
            lists (e.g.  [[start, stop, node1, node2], ..]) event_keys must map
            to the corresponding indices (so {'start':0, 'stop': 1, ...}).
            Equivalently, if events is a list of dict then event_keys must map
            to the corresponding keys.
        :param source: either the path to a text file (e.g. csv) or a python
            dict holding all the events.
        :param params:
            :param key_mapping: optional dict to map custom keys of the contact
            sequence to the required keys, that are 'start', 'stop', 'node1'
            and 'node2'
        :return:
    """
    def __init__(
            self, nodes=None, node_import=None, events=None, event_keys=None,
            source=None, **params
            ):
        # _Graph.__init__(self)
        self.t_start, self.t_stop = None, None
        self._nodes_start, self._nodes_stop = None, None
        self.nodes_start, self.nodes_stop = None, None
        self.is_static = False
        self.has_dynamic_nodes = False
        # make type specific imports
        if nodes is not None:
            self.nodes = nodes
            self.o_ids = [node.uid for node in self.nodes]
            n = len(self.o_ids)
            self.nodes_start = np.array([node.start for node in self.nodes])
            self.nodes_end = np.array([node.stop for node in self.nodes])
            if node_import is not None:
                for attr in node_import:
                    setattr(
                            self, attr,
                            [getattr(node, attr, None) for node in self.nodes]
                            )
            # ToDo: deal with the nodes contacts
            if events is not None:
                events.sort(key=lambda x:x[event_keys['start']])
                self._node1s = [
                    an_event.get(event_keys['node1']) for an_event in events
                    ]
                #self.node1s = [
                #    Node.get_node(node)._id for node in self._node1s
                #    ]
                self._node2s = [
                    an_event.get(event_keys['node2']) for an_event in events
                    ]
                #self.node2s = [
                #    Node.get_node(node)._id for node in self._node2s
                #    ]
                # NOTE: other stage had get for self.starts and self.stops
                self.starts = np.array([
                    an_event.pop(event_keys['start']) for an_event in events
                    ])
                self.stops = np.array([
                    an_event.pop(event_keys['stop']) for an_event in events
                    ])
                self.event_params = events
        elif events is not None:
            # ToDo: Handle the import from events only case
            pass
        else:
            if isinstance(source, str):  # it is a file
                self._load_from_file(source, **params)
            elif isinstance(source, dict):
                # if we don't need to pass on the params, better don't
                self._load_from_dict(source, **params.get('key_mapping', {}))
            else:  # source must be an EventQueue then
                # todo: copy events that were passed by arguments
                raise NotImplementedError('Provide either the path to a file' +
                        ' or a dictionary as source')
            # ToDo: self.event_params is not filled
            self.event_params = []
            self.o_ids = list(np.union1d(self._node1s, self._node2s))
        # remap the node ids
        n = len(self.o_ids)
        mapper = {val: key for key, val in enumerate(self.o_ids)}
        get_element = lambda k: mapper.get(k)
        v_get_id = np.vectorize(get_element)

        # self.node1/2s is the list of 'usable' node ids.
        # self._node1/2s are the original node ids
        self.node1s = v_get_id(self._node1s)
        self.node2s = v_get_id(self._node2s)

        # If nodes_start and nodes_end are specified in params, overwrite it
        if 'nodes_start' in params:
            self._nodes_start = params.get('nodes_start')
        if 'nodes_end' in params:
            self._nodes_end = params.get('nodes_end')

        if self._nodes_start not in [{}, None]:
            self.nodes_start = {
                    v_get_id(node): self._nodes_start[node]
                    for node in self._nodes_start
                    }
        if self._nodes_end not in [{}, None]:
            self.nodes_end = {
                    v_get_id(node): self._nodes_end[node]
                    for node in self._nodes_end
                    }

        self.t_start = params.get(
                't_start',
                np.min(self.starts) if self.t_start is None else self.t_start
                )
        self.t_stop = params.get(
                't_stop',
                np.max(self.stops) if self.t_stop is None else self.t_stop
                )
        _Graph.__init__(self, n=n)


        # transform nodes_start and nodes_end into np.arrays if needed
        if not (self.nodes_start is None or self.nodes_end is None):
            if len(self.nodes_start) != n or len(self.nodes_end) != n:
                InvalidArgumentError(
                        'The <nodes_start> and <nodes_end> arguments have to '\
                                'be of the same length as the total number of '\
                                'unique nodes (n=%s)', str(self.n))
            self._transform_nodes_start_end(v_get_id)
            self.has_dynamic_nodes = True
        else:
            # set to default values
            self.nodes_start = np.repeat(self.t_start, n)
            self.nodes_end = np.repeat(self.t_stop, n)

        _Graph.__init__(self, n=n)

    def _transform_nodes_start_end(self, vectorize_func):
        """
        Transform self.nodes_start and self.nodes_end into numpy arrays in
        which the value at position i is the value of individual i, defined by
        the vectorize function. If self.nodes_start and self.nodes_end are
        already numpy arrays, remap positions

        Arguments:

        :param vectorize_func: Maps an ID of general type (e.g. String) to an
            integer in the space {0...n}, where n is the number of nodes
        :return:
        """
        if isinstance(
                self.nodes_start, dict
                ) and isinstance(
                        self.nodes_end, dict
                        ):
            nodes_start = copy(self.nodes_start)
            nodes_end = copy(self.nodes_end)
            self.nodes_start = np.zeros(len(nodes_start))
            self.nodes_end = np.zeros(len(nodes_end))
            for node_name, val in nodes_start.items():
                self.nodes_start[vectorize_func(node_name)] = val
            for node_name, val in nodes_end.items():
                self.nodes_end[vectorize_func(node_name)] = val

        elif isinstance(
                self.nodes_start, np.ndarray
                ) and isinstance(
                        self.nodes_end, np.ndarray
                        ):
            # in this case self.o_ids has to be a list of integers with all
            # integer values up to n that map to positions in nodes_start and
            # nodes_end. Simply re-map positions...
            self.nodes_start = self.nodes_start[vectorize_func(self.o_ids)]
            self.nodes_end = self.nodes_end[vectorize_func(self.o_ids)]
        else:
            InvalidArgumentError(
                    'The arguments <nodes_start> and <nodes_end> have to be '\
                            'of type dict or numpy.array'
                            )


    def _load_from_file(self, source, **params):
        """

        Parameters:
        -----------
        :param source:
        :param params: Several arguments depending on the type of the source
            argument.
            Mandatory:
                start_tag: Name of the column containing the
                    start time/frame tag.
                stop_tag: Name of the column containing the
                    stop time/frame tag.
                node1_tag: Name of the column with the first
                    participant
                node2_tag: Name of the column with the second
                    participant
                delimiter: The string delimiting the columns

                    (default: 'TAB'). You can use the actual
                    string e.g. '\\t' for 'TAB" or choose from:
                    ('TAB', 'space')
                string_values: A list of all the columns containing
                    strings as values. All the
                    other columns will be converted to floats.
            Optional:
                directed: True/False whether the interactions
                        are directed or not.

            If source is a filename
                Optional:
                    permitted_values: Default = {}
                        If 'start_tag' is given, only interactions with

                        a start_tag bigger than this value
                            will be considered.

            If source is an EventQueue:
                Optional:
                    t_start: The time/frame at which the scenario
                        starts.
                    t_end: The time/frame at which the scenario
                        ends.
                    nodes: List of considered nodes.
                        Note: all other nodes are ignored as well as
                        interactions with them.
        :return:
        """

        try:
            with open(source, 'r') as f:
                pass
        except IOError:
            InvalidArgumentError('The file specified does not exist')
        for arg in [
                'start_tag', 'stop_tag', 'node1_tag', 'node2_tag',
                'string_values', 'delimiter'
                ]:
            try:
                if arg == 'delimiter':
                    if params[arg] == 'TAB':
                        delim = '\t'
                    elif params[arg] == 'space':
                        delim = ' '
                    else:
                        delim = params[arg]
                    self._file_delimiter = delim
                setattr(self, arg, params[arg])
            except KeyError:
                raise InvalidArgumentError(
                        'The mandatory argument <%s> was not given.' % arg
                        )
        # each event will have this structure
        self._event_structure = [
                self.start_tag, self.stop_tag, self.node1_tag, self.node2_tag
                ]
        with open(source, 'r') as f:
            self._file_header = f.readline().rstrip().replace(
                    '#', ''
                    ).split(self._file_delimiter)
            data = np.genfromtxt(
                f,
                delimiter=self._file_delimiter,
                unpack=False,
                autostrip=True,
                comments='#',
                names=', '.join(self._file_header),
                # try to determine type independently
                dtype=None,
                usecols=tuple(
                    filter(
                        lambda x: x in self._event_structure, self._file_header
                        )
                    )
                )
            self.starts = data[self.start_tag]
            self.stops = data[self.stop_tag]
            self._node1s = data[self.node1_tag]
            self._node2s = data[self.node2_tag]
            # ToDo: Read nodes_start and nodes_end as optional arguments from a
            # text file

    def _load_from_dict(self, source, **params):
        """
        Create a temporal graph from a python dictionary. The following keys
            are mandatory:

            - start: a list/array of start times for each event
            - stop: a list/array of end times for each event
            - node1: a list/array of IDs for interaction partner 1 for each
                event
            - node2: a list/array of IDs of interaction partner 2 for each
                event

            Optional:
            - t_start: Start of simulation (default: start time of earliest
                event)
            - t_end: End of simulation (default: end time of last event)
            - nodes_start/ nodes_end: Time of start and end of the life-span of
                every node in the network. Allowed
            types:
                dictionary: Keys correspond to the node IDs, Values to the time
                    of start/end
                array: Value at index i corresponds to the time of start/end of
                    individual i

        Parameters:
        -----------

        :param source:
        :return:
        """

        _start = params.get('start', 'start')
        _stop = params.get('stop', 'stop')
        _node1 = params.get('node1', 'node1')
        _node2 = params.get('node2', 'node2')
        _tstart = params.get('t_start', 't_start')
        _tstop = params.get('t_stop', 't_stop')
        _nstart = params.get('first', 'first')
        _nend = params.get('last', 'last')
        # mandatory arguments:
        try:
            self.starts = np.array(source[_start], dtype=np.float64)
            self.stops = np.array(source[_stop], dtype=np.float64)
            self._node1s = np.array(source[_node1])
            self._node2s = np.array(source[_node2])
        except KeyError:
            raise InvalidArgumentError(
                    'Loading the temporal graph from a dict failed.\n Here is'
                    'how to do this porperly:\n\n%s' % (
                        self._load_from_dict.__doc__
                        )
                    )

        # Optional part
        self.t_start = np.array(source.pop(_tstart, np.min(self.starts)))
        self.t_stop = np.array(source.pop(_tstop, np.max(self.stops)))
        # get node start and stop
        self._nodes_start = {
                _node: start
                for _node, start in params.get('nodes_start', {}).items()
                }
        self._nodes_end = {
                _node: end
                for _node, end in params.get('nodes_end', {}).items()
                }


        # If further data for the nodes is present, pass them to self.params
        self.host_params = source.get('host_params', {})


    # ToDo: will be replaced by _load_from_dict
    def _copy_events(self, **params):
        """ copy events informations from existing arrays given as keyword
            arguments.

        Parameters
        ----------

        starts: float array
            starting times of the meetings
        stops: float array
            stoping times of the meetings
        node1s: int array
            IDs of the first mice for each meetings
        node2s: int array
            IDs of the second mice for each meetings
        """

        self.starts = np.array(params['starts'], dtype=np.float64)
        self.stops = np.array(params['stops'], dtype=np.float64)
        self.node1s = np.array(params['node1s'], dtype=np.int64)
        self.node2s = np.array(params['node2s'], dtype=np.int64)


class Graph(_Graph):
    def __init__(self, n=None, method='stub', **distribution):
        """
            Possible arguments for the distribution are:
            - network_type: specify the type of network that should be
                constructed (THIS IS MANDATORY).
                It can either be the name of a distribution or of a certain
                    network type.

            ['uniform', 'full', 'l_partition', 'poisson', 'normal', 'binomial',
                'exponential', 'geometric', 'gamma', 'power', 'weibull']

            For specific parameters of the distributions, see:
                http://docs.scipy.org/doc/numpy/reference/routines.random.html

            - method: The probabilistic framework after which the network will
                be constructed.
            - distribution specific arguments. Check out the description of the
                specific numpy function. Or just give the argument network_type
                and look at what the error tells you.
                
           See self._create_graph for more information
        """
        _Graph.__init__(self)
        self.is_static = True
        self.has_dynamic_nodes = False
        self._rewiring_attempts = 100000
        self._stub_attempts = 100000
        self.permitted_types = allowed_dists + [
                "l_partition", 'full', 'uniform'
                ]
        self.is_directed = False
        # to do: pass usefull info in here.
        self._info = {}
        # for now only undirected networks
        if n is not None:
            self.n = n
            if method in ['proba', 'stub']:
                self.method = method
            else:
                raise ValueError(
                        method + ' is not a permitted method! Chose either \
                                "proba" or "stub"'
                                )
            try:
                self.nw_name = distribution.pop('network_type')
                empty_graph = False
            except KeyError:
                self.nn = []
                self._convert_to_array()
                empty_graph = True
                # create an empty graph if network_type is not given
            if not empty_graph:
                if self.nw_name not in self.permitted_types:
                    raise ValueError(
                        "The specified network type \"%s\" is not permitted. \
                        Please chose from " % self.nw_name + '[' + ', '.join(
                            self.permitted_types) + ']'
                        )
                self.distribution = distribution
                self._create_graph(**self.distribution)

            # create the o_ids attribtue
            self.o_ids = range(self.n)

    def _create_graph(self, **distribution):
        """
            Creates an explicit graph.
            The degrees are drawn for the wanted distribution.
            The graph is then created using either the stub-algorithm
                or with a probabilistic method (by the use of self._make_graph.
            The stub-algorithm keeps the exact degree that was drawn from
                the distribution, the probabilistic method on the other side
                produces a realization where the drawn degree is only the
                expected value. 
            Arguments are:
                method: either 'stub' or 'proba', default = 'stub'
                    - determines which algorithm is used
                    Note that the stub algorithm will have an increasingly hard
                        time for denser networks and will eventually fail if
                        the average degree is almost n-1.
            **distribution: several arguments.
                network_type: gives the name of the wanted distribution
        """
        if self.nw_name == 'l_partition':
            try:
                self.l_partition_network(**distribution)
            except TypeError, msg:
                print 'OH, something went wrong! \n Here, have a description \
                        of the distribution:'
                print self.l_partition_network.__doc__
                raise TypeError(msg)
            return 0
            #is needed since it should not run through the rest in this case.
        elif self.nw_name == 'full':
            self.fully_connected(**distribution)
            return 0
        elif self.nw_name == 'uniform':
            try:
                degrees = [distribution['degree'] for _ in xrange(self.n)]
            except KeyError, msg:
                print 'You need to provide the argument degree if you \
                        want to create a network with uniform degree!'
                raise KeyError(msg)
        else:
            print distribution
            try:
                degrees = []
                for _ in xrange(self.n):
                    degrees.append(Distribution[self.nw_name](**distribution))
            except TypeError, msg:
                print 'OH, something went wrong! \n Here, have a description '\
                        'of the distribution:'
                print Distribution[self.nw_name].__doc__
                raise TypeError(msg)
        #except KeyError:
        #    try:
        #        degrees = []
        #        for _ in xrange(self.n):
        #            degrees.append(Distribution[self.nw_name](**distribution))
        #    except TypeError, msg:
        #        print 'OH, something went wrong! \n Here, have a description \
        #            of the distribution:'
        #        print Distribution[self.nw_name].__doc__
        #        raise TypeError(msg)
        self._make_graph(degrees)
        return 0

    class ConstructionError(Exception):
        pass

    # to do: look at the stub algorithm, something is not working correctly
    def _make_graph(self, degrees, stub_tries=100):
        """
            This function constructs the actual network.
            Depending on the construction method that has been chosen the

                network is either generated by a probabilistic method
                (if self.method = 'proba') or by a edge wiring method
                (self.method = 'stubs').
        """
        if self.method == 'proba':
            two_k = float(sum(degrees))
            self.nn = [[] for _ in xrange(self.n)]
            for i in xrange(len(degrees) - 1):
                for j in xrange(i + 1, len(degrees)):
                    if min(
                            degrees[i] * degrees[j] / two_k, 1
                            ) > random.sample():
                        self.nn[i].append(j)
                        self.nn[j].append(i)
            self._convert_to_array()
            self.degrees = map(lambda x: len(x), self.nn)
            return 0
        else:
            if sum(degrees) > self.n * (self.n - 1):
                print 'WARNING: according to the chose degree sequence the \
                        network should be more than fully connected'
                print '\t the degree sequence is corrected such that each node \
                        has at most n-1 neighbours.'
            self.degrees = map(
                    lambda x: min(self.n - 1, int(round(x))), degrees
                    )
            # check if the sum of the degrees is odd
            s_deg = sum(self.degrees)
            # if odd, correct (+1) as otherwise the stub method cannot work
            if s_deg/2 != s_deg/2.:  
                self.degrees[0] += 1
            stubs = []
            self.nn = [[] for _ in xrange(self.n)]
            for node in xrange(self.n):
                stubs.extend([node for _ in xrange(self.degrees[node])])
            # do self._stub_attempts attempts to get the connections right
            while stub_tries:
                length = len(stubs)
                length_queue = Queue(maxsize=self._stub_attempts)
                for _ in xrange(self._stub_attempts):
                    length_queue.put_nowait(0)
                while length:
                    if length/2 != length/2.:
                        print length, 'on'
                    try:
                        n_1 = self._get_rand_element(stubs)
                    except ValueError:
                        break
                    try:
                        n_2 = self._get_rand_element(stubs)
                    except ValueError:
                        stubs.append(n_1)
                        break
                    if n_1 != n_2:
                        if n_2 not in self.nn[n_1]:
                            self.nn[n_1].append(n_2)
                            self.nn[n_2].append(n_1)
                        else:
                            stubs.append(n_2)
                            stubs.append(n_1)
                    else:
                        stubs.append(n_2)
                        stubs.append(n_1)
                    # need only to be recomputed if stubs.append is not executed 
                    length = len(stubs)  
                    # if after self._stub_attempts the stubs has the same
                    if length_queue.get_nowait() - length == 0:  
                        # length, give up
                        break
                    length_queue.put_nowait(length)
                #here the remaining elements from the stubs list
                #cannot be matched.
                counter = copy(self._rewiring_attempts)
                if len(stubs):
                    print 'trying to fix'
                    print len(stubs)
                while len(stubs) and counter:
                    n_1 = self._get_rand_element(stubs)  # take a stubs
                    # get the list of potential neighbours for this stub
                    pot_neigh_n_1 = filter(
                            lambda x: x not in self.nn[n_1], range(self.n)
                            )
                    try:
                        rand_node_1 = self._get_rand_element(pot_neigh_n_1)
                        # get a neighbour of rand_node_1
                        rand_node_2 = self.nn[
                                rand_node_1
                                ][random.randint(0, len(self.nn[rand_node_1]))]
                    except ValueError:
                        stubs.append(n_1)
                        continue
                    self.nn[n_1].append(rand_node_1)  
                    # connect the stub to the random node
                    self.nn[rand_node_1].append(n_1)  # "
                    self.nn[rand_node_1].remove(rand_node_2)  #
                    self.nn[rand_node_2].remove(rand_node_1)
                    # rand_node_2 has a free stub now try to match it with
                    # another stub
                    candidates = filter(
                            lambda x: x not in self.nn[rand_node_2], stubs
                            )
                    if candidates:
                        new_mate = self._get_rand_element(candidates)
                        # remove the new_mate from the stubs
                        stubs.pop(stubs.index(new_mate))  
                        self.nn[rand_node_2].append(new_mate)
                        self.nn[new_mate].append(rand_node_2)
                    else:
                        # nothing to match, put the stub of rand_node_2 into
                        # the list
                        stubs.append(rand_node_2)
                    print len(stubs), counter
                    counter -= 1
                if not len(stubs):  # it is constructed
                    self._convert_to_array()
                    return 0
                print 'new_try'
                stub_tries -= 1
            raise self.ConstructionError("""The Graph could not be constructed.
            You can try to launch the function again, but
            consider that you probably chose a graph that is
            too dense which makes it hard for the stubs algorithm
            to find a coherent list of links. In case it fails several
            times consider using the "proba" method."""
            )

    @property
    def info(self):
        """
            Returns the set of parameters used to create the network.
        :return: Dictionary with the network type and the set of parameters
            used to generate the network
        """
        return {
            'n': self.n,
            'constuction_method': self.method,
            'type': self.nw_name,
            'type_specific': self.distribution

        }

    def fully_connected(self, **kwargs):
        """
        Create a fully connected network
        :param kwargs:
        :return:
        """
        self.nn = []
        for node in xrange(self.n):
            nns = range(self.n)
            nns.remove(node)
            self.nn.append(nns)
        self.degrees = [self.n - 1 for _ in xrange(self.n)]
        self._convert_to_array()

    def l_partition_network(self, l, **kwargs):
        """
        Create an l_partition graph.
        :param l: number of partitions
                Note: if n is not a multiple of l,
                    some partitions will be smaller/bigger
                    than others.
        :param kwargs: Set of possible arguments to define an l-partition
            network:
                Possible arguments are:
                    - avg_degree: give the average degree
                - density_ratio: the edge density ratio between
                    inside and in-between partitions. A ratio
                    of 2 means inside is twice as dense as outside
                - p_in: the connection probability for a pair
                    inside a partition.
                - p_out: the connection probability for a pair
                    of nodes between partitions.
        """
        if 'p_in' in kwargs and 'p_out' in kwargs:
            p_in = kwargs['p_in']
            p_out = kwargs['p_out']
        elif 'avg_degree' in kwargs and 'density_ratio' in kwargs:
            p_in = kwargs['avg_degree'] / float(
                    (
                        self.n / float(l) - 1
                        ) * (
                            1 + 1 / float(kwargs['density_ratio'])
                            )
                    )
            p_out = kwargs['avg_degree'] / float(
                    (l - 1) * self.n / float(l) * (1 + kwargs['density_ratio'])
                    )
            if max(p_in, p_out) > 1:
                raise ValueError(
                    """The desire ratio along with the given
                    average degree is not a possible realization.
                    For this to work, parts of the network would need to
                    be more than fully connected."""
                )
        else:
            raise TypeError(
                    'Invalid arguments for the l_partition_network function'
                    )
        nodes = xrange(self.n)
        partition_size = int(self.n // float(l))
        part_sizes = [partition_size for i in xrange(l)]
        rest = self.n % l
        while rest:
            part_sizes[rest - 1] += 1
            rest -= 1
        partitions = [0]
        for a_size in part_sizes:
            partitions.append(partitions[-1] + a_size)
        partitions = zip(partitions[0:-1], partitions[1:])
        self.nn = [[] for i in nodes]
        for i in xrange(len(partitions)):
            a_partition = partitions[i]
            for node_1 in xrange(a_partition[0], a_partition[1] - 1):
                #inside partition node connections
                for int_node_2 in xrange(node_1 + 1, a_partition[1]):
                    if random.sample() <= p_in:
                        self.nn[node_1].append(int_node_2)
                        self.nn[int_node_2].append(node_1)
                #external links
                if l - i - 1:
                    for ext_node_2 in xrange(a_partition[1], self.n):
                        if random.sample() <= p_out:
                            self.nn[node_1].append(ext_node_2)
                            self.nn[ext_node_2].append(node_1)
        self.degrees = [len(nn) for nn in self.nn]
        self._convert_to_array()
        return 0

    def _convert_to_array(self):
        for i in xrange(len(self.nn)):
            self.nn[i] = array(self.nn[i])

    @staticmethod
    def _get_rand_element(a_list):
        """
            Returns and removes a random element from a list.
        """
        rand_index = random.randint(0, len(a_list))
        #swap element with last one (to use pop on it)
        a_list[rand_index], a_list[-1] = a_list[-1], a_list[rand_index]
        return a_list.pop()

    def _get_nn(self, edge_list):
        """
            Construct a list of nearest neighbors (nn) for each node.
            
        """
        self.nn = [[] for _ in xrange(self.n)]
        for edge in edge_list:
            self.nn[edge[0]].append(edge[1])
            self.nn[edge[1]].append(edge[0])
        self._convert_to_array()

    def new_realisation(self, **distribution):
        """
            Creates a new realisation of the chosen graph type.
            If no arguments are specified previously defined parameters are
                used.
        """
        if distribution:
            for key in distribution:
                self.distribution[key] = distribution[key]
        self._create_graph(**self.distribution)
        return 0

    def _get_edgelist(self, ):
        """
        Used to compute the edge list for the export to file.
        """
        edgelist = []
        for node in xrange(self.n):
            for nn in self.nn[node]:
                edge = ' '.join([str(node), str(nn)]) + '\n'
                if edge not in edgelist and ' '.join(
                        [str(nn), str(node)]
                        ) + '\n' not in edgelist:
                    edgelist.append(edge)
        return edgelist

    def export_graph(self, filename, fileformat=None):
        """
        This function exports the graph as it is to a file.
        Arguments:
            - filename: The name of the file to hold the graph.
                eg.: '/home/user/FancyGraphs/my_fancygraph.txt'
                Note: If the filetype is given, the fileformat argument
                    does not need to be set.
            - fileformat: Default is None. If no format is set, the format is 
                inferred from the filename.
                Permitted formats are:
                    *.txt/*.edges: Edgelist. Each column is a pair of integers
                        separated by whitespace
                    *.ncol: Large Graph Layout Format. Each column is a pair of
                        node names and an optional weight(pos/neg integer),
                        each separated by whitespace
                    *.lgl: Large Graph Layout (alternative format). Adjacency
                        style type of list (so a list of lists) indicating an
                        source node with # followed by a list of target nodes
                        and optional weights.
                        Eg.
                        #node1
                        node2 -1
                        node4 
                        #node2
                        node3 10
        """
        ok_formats = ['txt', 'ncol', 'lgl', 'edges']
        if '.' in filename:
            [filename, fileformat] = filename.split('.')
        if fileformat not in ok_formats:
            raise IOError(
                    'The chosen format <%s> is not permitted. See docstring \
                            for info.' % fileformat
                            )
        with open(filename + '.' + fileformat, 'w') as the_file:
            if fileformat == 'lgl':
                for node in xrange(self.n):
                    the_file.write('#' + str(node) + '\n')
                    for nn in self.nn[node]:
                        #no weights for now.
                        the_file.write(str(nn) + '\n')
            elif fileformat == 'ncol':
                #at the moment weighted edges are not supported.
                edgelist = self._get_edgelist()
                the_file.write(''.join(edgelist))
            else:
                edgelist = self._get_edgelist()
                the_file.write(''.join(edgelist))
        return 0

    def import_graph(self, filename):
        """
        This function import a graph from a file.
        """
        ok_formats = ['txt', 'ncol', 'lgl', 'edges']
        extension = filename.split('.')[-1]
        if extension not in ok_formats:
            raise IOError('The file is of an unsupported format. Please \
                    provide a file with one of the following extensions:\n \
                    %s' % ', '.join(ok_formats))
        with open(filename, 'r') as the_file:
            lines = the_file.readlines()
            if extension in ['txt', 'edges', 'ncol']:
                edge_list = []
                for line in lines:
                    edge_list.append(
                            tuple(map(lambda x: int(x), line.split(' ')[:2]))
                            )
                    #only add the first two, other columns are weight etc.
                    #NOTE: The weights are not handled yet.
                self._get_nn(edge_list=edge_list)
            else:
                self.nn = []
                for line in lines:
                    if line[0] == '#':
                        key = int(line[1:])
                        #issue: catch exceptions
                        self.nn[key] = []
                    else:
                        elements = line.split(' ')
                        try:
                            self.nn[key].append(int(elements[0]))
                        except NameError:
                            raise NameError(
                            'The content needs to start with a #nodeID to \
                                    specify to whom the following nearest \
                                    neighbours belong'
                                    )
                        #issue: consider possible edge weights
                        #issue: catch exceptions
                self._convert_to_array()
            self.n = len(self.nn)
        return 0

    @property
    def as_igraph(self):
        """
        Method to export the Graph object to an igraph.Graph object
        :return: igraph.Graph
        """
        from igraph import Graph as iGraph
        as_igraph = iGraph(directed=False)
        as_igraph.add_vertices(range(self.n))
        for i in xrange(len(self.nn)):
            for neighbour in self.nn[i]:
                as_igraph.add_edges([(i, neighbour)])
        as_igraph.simplify()
        return as_igraph

    @property
    def as_networkx(self):
        """
        Method to export the Graph object ot a networkx.Graph object
        :return: networkx.Graph
        """
        from networkx import Graph as nxGraph
        as_nxgraph = nxGraph()
        as_nxgraph.add_nodes_from(xrange(self.n))
        for i in xrange(len(self.nn)):
            for neighbour in self.nn[i]:
                as_nxgraph.add_edge(i, neighbour)
        return as_nxgraph
