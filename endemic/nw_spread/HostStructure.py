__author__ = 'Jonas I Liechti'
import numpy as np


class ContactStructure():
    class HostOrderError(Exception):
        pass

    class IncompatibleSusceptibilityError(Exception):
        pass

    class UniqueIDError(Exception):
        pass

    def __init__(self, from_object, susceptible=1, is_static=True):
        self.is_static = is_static  # whether this is a static or a dynamic (time explicit network)
        suscept_default = 1  #if any susceptibility information is missing, it will be completed with this value.
        if isinstance(from_object, list):  # assume it is a list of Host objects:
            len_hosts = len(from_object)
            self._hosts = []
            self.nn = [[] for _ in xrange(len_hosts)]
            self._susceptible = [{} for _ in xrange(len_hosts)]
            for a_host in from_object:
                self._hosts.append(a_host.ID)
            self._hosts.sort()
            self.n = len(self._hosts)  #gives the number of hosts
            self._check_integrity()
            id_map = {}
            for val in self._hosts:
                id_map[val] = self._hosts.index(val)
            for a_host in from_object:
                its_id = id_map[a_host.ID]
                self.nn[its_id] = a_host.neighbours
                its_default = suscept_default
                susceptibility = a_host.susceptible
                if 'Default' in susceptibility:
                    its_default = susceptibility.pop('Default')
                self._susceptible[its_id]['Default'] = its_default
                for strain_name in susceptibility.keys():
                    self._susceptible[its_id][strain_name] = susceptibility[strain_name]
                    #self.susceptible[its_id] = a_host.susceptible
            _hosts = id_map.values()
            _hosts.sort()
            self._hosts = _hosts
        elif True:  # it is a _Graph object # to do: set condition on being of _Graph class (with super?)
            self.is_static = from_object.is_static
            self.graph_info = from_object.info
            self._hosts = range(from_object.n)
            self.n = from_object.n
            if self.is_static:
                self.nn = from_object.nn
            else:
                self.nn = None
            if type(susceptible) is dict:
                if 'Default' not in susceptible:
                    susceptible['Default'] = suscept_default
                self._susceptible = [susceptible for _ in xrange(from_object.n)]
            elif type(susceptible) is float or type(susceptible) is int:
                self._susceptible = [{'Default': susceptible} for _ in xrange(from_object.n)]
            elif susceptible is list:
                self._susceptible = []
                for a_state in susceptible:
                    self._susceptible.append(a_state)
            else:
                raise self.IncompatibleSusceptibilityError(
                    """Susceptible must either be a list a dict or a float/int. Please refer to the class description\
                     for more details.""")
        else:
            raise AttributeError(
                'Neither a graph nor a lists of hosts is provided.\n\n%s' % self.__doc__
            )
        self.graph_info = {}
        # this is only filled up in the Scenario class in the Spreading module.
        self.susceptible = [[] for _ in xrange(len(self._susceptible))]

    @property
    def info(self):
        # ToDo: would make sense to define this.
        return self.graph_info

    def _check_integrity(self):
        """
        Check the integrity of the ensemble of hosts
        """
        if len(list(set(self._hosts))) != len(self._hosts):
            raise self.UniqueIDError(
                'Not all hosts have distinct IDs. Please provide a set of host with all unique name parameters'
            )
        else:
            pass
        return 0


class ContactNetwork(ContactStructure):
    def __init__(self, hosts=None, graph=None, susceptible=1):
        """
        The contact_structure class defines the appropriate Network of _hosts on which
            pathogens will spread.
        It either takes a graph from nw_construct package or a list of hosts as
            an argument. Note that if a list of _hosts is provided, they need to
            have the neighbours argument filled, otherwise no network is constructed.
            
        :type susceptible: int, dict, list, float
        Arguments:
            - graph: an object from the graph class defined in the nw_construct package.
            
            - hosts: a list of Host objects.
                Note: If hosts is provided the argument 'susceptible' is ignored.
            
            Note: One of the two arguments must be provided. If both are, the graph argument
                is ignored.

            - susceptible: Determines the susceptibility of hosts.
                Can be the status of a single node. If susceptible is a float [0,1] then for all strains for all
                nodes, this state is taken. If it is a single dict, then the list is taken for all nodes, i.e. all nodes
                have the same susceptibility for the strains. If it is a list of dicts, then each node gets its own
                susceptibility state, i.e. the index in the list determines the name of the node.
                Eg. A single dict: susceptible={'wild_type':1, 'resistant_1': 0, 'Default': 1}
                    which reads: susceptible to the 'wild_type' strain, resistant to 'resistant_1' and susceptible to
                    all other strains ('Default').
        """
        ContactStructure.__init__(self, from_object=graph if graph else hosts, is_static=True)

    def get_events(self, node_id):
        return self.nn[node_id]

    def update_topology(self, graph):
        """
        This method takes a new topology and updates the contact network accordingly.
        Note: the parameter graph must be a Graph from nw_construct module that shares the same properties
        as the previousely defined contact_structure, with the exception of the node linking.

        Usage e.g.:
            my_graph = nw_construct.Graph(...)
            my_contact_nw = ContactNetwork(graph=my_graph)
            #do simulations etc.
            #we want to update the topology:
            my_graph.new_realization()  #change my_graph
            my_contact_nw.update_topology(graph=my_graph)  #update my_coontact_nw according to new topology of my_graph

        :param graph:
        :type graph: nw_construct.Graph
        :return:
        """
        self.nn = graph.nn
        return None


class ContactSequence(ContactStructure):
    def __init__(self, temporal_graph=None, **params):
        self.starts = temporal_graph.starts
        self.stops = temporal_graph.stops
        self.node1s = temporal_graph.node1s
        self.node2s = temporal_graph.node2s
        self.t_start = params.get('t_start', np.min(self.starts))
        self.t_stop = params.get('t_stop', np.max(self.stops))
        ContactStructure.__init__(self, from_object=temporal_graph, is_static=False)

    def get_events(self, node_id, start_time, delta_t):
        """
        Returns a view of the start times and stop times as well as the involved nodes of all event vor a given
        node within a time range (start_time, start_time + delta_t)
        :param node_id:
        :param start_time:
        :param delta_t:
        :return:
        """
        stop_time = start_time + delta_t
        the_filter = np.logical_and(
            np.logical_and(
                start_time <= self.stops,
                stop_time > self.starts
            ), np.logical_or(
                node_id == self.node1s,
                node_id == self.node2s
            )
        )
        nn1 = self.node1s.view()[the_filter]
        nn2 = self.node2s.view()[the_filter]

        nn = np.where(
            node_id != nn1,
            nn1,
            nn2
        )
        return nn, self.starts.view()[the_filter], self.stops.view()[the_filter]


class Host():
    def __init__(self, an_id, neighbours, susceptible):
        """
                This class defines a single host.

        Arguments:
            - name: int, unique for each host in a contact_structure
            - neighbours:
            - susceptible:
        :param an_id: unique id for each host in a contact_structure
        :param neighbours: A numpy array of either host _id's indicating all the neighbours of the host or tuples
        :param susceptible: A float [0,1] or dict. If it's a dict, the key is the name of a strain and the corresponding
                value [0,1] indicates whether the host is susceptible or not. Additionally the key 'Default' can be
                given. It will associate to all missing strains the specified value.
                If a float is given, the value will be taken for all strains. This is equivalent to just give a
                    dict with a 'Default' value.
        :return:
        """
        self._id = an_id
        self.susceptible = susceptible
        self.neighbours = neighbours