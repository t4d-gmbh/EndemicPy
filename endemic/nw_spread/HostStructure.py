__author__ = 'Jonas I Liechti'
import numpy as np
from Queue import Queue

class ContactStructure():
    class HostOrderError(Exception):
        pass

    class IncompatibleSusceptibilityError(Exception):
        pass

    class UniqueIDError(Exception):
        pass

    class ImplementationMissingError(Exception):
        pass

    def __init__(self, from_object, susceptible=1, is_static=True, has_dynamic_nodes=False):
        self.is_static = is_static  # whether this is a static or a dynamic (time explicit network)
        self.has_dynamic_nodes = has_dynamic_nodes  # in the case of a temporal graph nodes can also be dynamic
        suscept_default = 1  # if any susceptibility information is missing, it will be completed with this value.
        self.o_ids = []  # match between the original ids of the nodes (values) and the newly generated ones (index)
        if isinstance(from_object, list):  # assume it is a list of Host objects:
            if self.is_static:
                len_hosts = len(from_object)
                self._hosts = []
                self.nn = [[] for _ in xrange(len_hosts)]
                self._susceptible = [{} for _ in xrange(len_hosts)]
                for a_host in from_object:
                    self._hosts.append(a_host.ID)
                self._hosts.sort()
                self.n = len(self._hosts)  # gives the number of hosts
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
                # create the original ID's list
                inv_id_map = {v: k for k, v in id_map.iteritems()}
                index_keys = inv_id_map.keys()
                index_keys.sort()
                self.o_ids = [inv_id_map[a_key] for a_key in index_keys]
            else:
                # ToDo: Create ContactStructure with a list of 'temporal' host objects.
                raise self.ImplementationMissingError("""Creating a ContactStructure with a list of 'temporal'
                Host objects is not implemented.""")

        elif True:  # it is a _Graph object # to do: set condition on being of _Graph class (with super?)
            self.is_static = from_object.is_static
            self.graph_info = from_object.info
            self.o_ids = from_object.o_ids
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
        self.nodes_start = temporal_graph.nodes_start
        self.nodes_end = temporal_graph.nodes_end
        ContactStructure.__init__(self, from_object=temporal_graph, is_static=False,
                                  has_dynamic_nodes=temporal_graph.has_dynamic_nodes)

    def get_events(self, node_id, start_time, delta_t):
        """
        Returns a view of the start times and stop times as well as the involved nodes of all event for a given
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
        nn1 = self.node1s[the_filter]
        nn2 = self.node2s[the_filter]

        nn = np.where(
            node_id != nn1,
            nn1,
            nn2
        )
        return nn, self.starts[the_filter], self.stops[the_filter]

    def get_nodes_by_lifetime(self, time_start, time_stop):
        """
        Returns all nodes that are "active" (as specified by TemporalGraph.nodes_start and TemporalGraph.nodes_end.

        :param time_start:
        :param time_stop:
        :return:
        """
        node_indices = np.logical_and(self.nodes_start < time_stop, self.nodes_end > time_start)

        return [i for i in xrange(len(node_indices)) if node_indices[i]]

    def get_temporally_connected_nodes(self, source_node, start_time, delta_t):
        """
        Finds all the nodes connected to source_node by a time-respecting path
        in a given time window
        
        Parameters
        ----------
        source_node : int
            the source node for the search        
        start_time : float 
            the starting time of the search
        delta_t: float
            stop_time = `start_time` + `delta_t`        
        
        Returns
        -------
        distances : list of ints
            distances from the source node for each node (-1 means that the node is unreachable)
        delays : list of floats
            time delays between each node and the source node
        
        Call
        ----
        distances, delays = get_temporally_connected_nodes(source_node, start_time, delta_t)
        
        11.08.2015, A.Bovet
        """
        # array holding the distance from the source node for all found nodes (-1 = not treated)
        distances = [-1 for _ in xrange(self.n)]    
        distances[source_node] = 0
    
        # arry holding the delay (temporal distance between the source node and the 
        # temporally connected nodes (-1 = not treated)
        delays = [-1 for _ in xrange(self.n)]    
        delays[source_node] = start_time    
        
        # queue containing the nodes whose neighbours need to be searched
        search_queue = Queue(maxsize = self.n)
        search_queue.put_nowait(source_node)
        
        while not search_queue.empty():
            node = search_queue.get_nowait()
            dist = distances[node]
            start = delays[node]
            # lookup time respecting neighbours         
            nn, nstarts, _ = self.get_events(node, start, start_time + delta_t - start)
            
            for neigh, neigh_start in zip(nn, nstarts):
                # if we haven't aready visited this node
                if distances[neigh] == -1:
                    distances[neigh] = dist + 1
                    delays[neigh] = neigh_start 
                    search_queue.put_nowait(neigh)
        
        delays = [delays[i]- start_time for i in xrange(len(delays))]
        return distances, delays
        
    def get_influence_set(self, node):
        """ returns the set of node indexes of the influence set of node "node"  """
        distances, _ = self.get_temporally_connected_nodes(node, self.t_start, self.t_stop - self.t_start)
        return [i for i, x in enumerate(distances) if x != -1]


# ToDo: Make temporal part of this class. (low priority)
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