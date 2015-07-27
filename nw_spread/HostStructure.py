class ContactNetwork():
    class HostOrderError(Exception):
        pass

    class IncompatibleSusceptibilityError(Exception):
        pass

    def __init__(self, hosts=None, graph=None, susceptible=1):
        """
        The contact_network class defines the appropriate Network of _hosts on which
            pathogens will spread.
        It either takes a graph from nw_construct package or a list of hosts as
            an argument. Note that if a list of _hosts is provided, they need to
            have the neighbours argument filled, otherwise no network is constructed.
            
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
        self.is_static = True  # set whether this is a static or a dynamic (time explicit network)
        suscept_default = 1  #if any susceptibility information is missing, it will be completed with this value.
        self.graph_info = {}
        if hosts:
            len_hosts = len(hosts)
            self._hosts = []
            self.nn = [[] for _ in xrange(len_hosts)]
            self._susceptible = [{} for _ in xrange(len_hosts)]
            for a_host in hosts:
                self._hosts.append(a_host.ID)
            self._hosts.sort()
            self.n = len(self._hosts)  #gives the number of hosts
            self._check_integrity()
            id_map = {}
            for val in self._hosts:
                id_map[val] = self._hosts.index(val)
            for a_host in hosts:
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
        elif graph:
            self.graph_info = graph._info
            self._hosts = range(graph.n)
            self.n = len(self._hosts)
            self.nn = graph.nn
            if type(susceptible) is dict:
                if 'Default' not in susceptible:
                    susceptible['Default'] = suscept_default
                self._susceptible = [susceptible for _ in xrange(graph.n)]
            elif type(susceptible) is float or type(susceptible) is int:
                self._susceptible = [{'Default': susceptible} for _ in xrange(graph.n)]
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
        self.susceptible = [[] for _ in xrange(len(self._susceptible))]
        # this is only filled up in the Scenario class in the Spreading module.

    def update_topology(self, graph):
        """
        This method takes a new topology and updates the contact network accordingly.
        Note: the parameter graph must be a Graph from nw_construct module that shares the same properties
        as the previousely defined contact_network, with the exception of the node linking.

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

    @property
    def info(self):
        return self.graph_info

    class UniqueIDError(Exception):
        pass

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


class Host():
    def __init__(self, ID, neighbours, susceptible):
        """
        This class defines a single host.
        
        Arguments:
            - name: int, unique for each host in a contact_network
            - neighbours: a numpy array of host ID's indicating all the neighbours of the host.
            - susceptible: A float [0,1] or dict. If it is a dict, the key is the name of a strain and the corresponding
                value [0,1] indicates whether the host is susceptible or not. Additionally the key 'Default' can be
                given. It will associate to all missing strains the specified value.
                If a float is given, the value will be taken for all strains. This is equivalent to just give a
                    dict with a 'Default' value.
        """
        self.ID = ID
        self.susceptible = susceptible
        self.neighbours = neighbours