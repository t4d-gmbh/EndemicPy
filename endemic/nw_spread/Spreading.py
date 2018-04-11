__author__ = 'Jonas I Liechti'
import random
from numpy import vectorize, array, where, absolute, zeros, count_nonzero
from numpy import random as nrand
from numpy import copy as n_copy
from numpy.ma import divide
from numpy import append as n_append
from numpy import sum as n_sum
from numpy import max as n_max
from Queue import Empty, PriorityQueue
from copy import copy, deepcopy
from collections import Counter
from RateDistribution import Distro


def _pick_rand_el(a_list):
    """
        Returns and removes a random element from a list.
    """
    rand_index = nrand.randint(0, len(a_list))
    # swap element with last one (to use pop on it)
    a_list[rand_index], a_list[-1] = a_list[-1], a_list[rand_index]
    # remove and return the element
    return a_list.pop()


def _get_rand_el(a_list):
    """
        Returns random element from a list without removing it.
    """
    return a_list[nrand.randint(0, len(a_list))]


class Scenario():
    def __init__(self, contact_structure, pathogen, treatment=None, **params):
        """
        Arguments:
            - contact_structure: An object of the contact_structure class (one of its subclasses), holding
                an ensemble of Host class objects.
                The spreading process will take place on this object.
            - pathogen: An object of the pathogen class, holding an
                ensemble of Strain class objects.
            - treatment: An object of the treatment class, holding an
                ensemble of drug class objects.
            - **params: Additional parameters for the simulation:
                - t_max: gives the maximal time after which the simulation stops.
                    Defaults=1000
                - t_burn: time window for the exclusive spread of the wild type.
                - dt: the time step at which the current status should be recorded.
                - default_susceptibility: Value [0,1] to be used for any missing susceptibility.
                    Default=1.
                - single_transmission: Boolean (default: False) indicating 
                    whether on contact between a carrying and a susceptible 
                    individual several transmission events should be possible or not.
                    Note: Traditional SIS, SIR or similar processes are not 
                    restrained to a single transmission per contact. 
                TODO: Is this even present?
                - ignore_dyn_nodes_in_log: If True, nodes that expand their lifespan are not set to -2 in self.log but
                    keep their state or they may undergo future state changes (e.g. recovery). This option only has an
                    effect if the contact_structure is a Temporal graph with specified node lifespans (arguments
                    nodes_start and nodes_end of the TemporalGraph instance)
        """
        # this will store self.current_view at various times
        # Note: this will be replaced with self.outcome so feel free to ignore this
        self.log = {}
        # this holds detailed information about what happened during a simulation
        self.simulation_log = {
            'phases': {},  # holds a dict with all the setup parameters at the starting time of a new phase.
            'mutations': [],  # holds all mutation events.
            'adjacency': {},  # the adjacency matrix of the contact network
            'modifications': {},  # keeps track of any change of the status (e.g. artificial infection event)
            'param_alternation': {}  # keeps track of any parameter alternations during the simulation
        }
        # this will hold processed output
        self.outcome = {}
        self.contact_structure = contact_structure
        # set the time:
        self.t = self.contact_structure.t_start if not self.contact_structure.is_static else 0
        self.pathogen = pathogen
        # init of optional and internal arguments
        self._default_susceptibility = 1  # by default hosts are susceptible
        self._default_drug_selection_factor = 1  # by default drugs do not increase mutation/selection rates
        # should a contact lead to multiple transmission events:
        self.single_transmission = params.get('single_transmission', False)
        # holds the number of infected individuals for each strain
        self._count_per_strains = array([0 for _ in xrange(self.pathogen.n)])
        # Note: used to contain the status of the host population overt time - will be redefined and is not in use.
        self._counts_over_time = zeros((1, self.pathogen.n))
        # set the default value for the time step used to regularly test if a quasi steady state is reached
        self._dt = params.get(
            'dt',  # issue: not sure yet about the default value
            5 / float(
                min(self.pathogen.trans_rates)
            ) if float(
                min(self.pathogen.trans_rates)
            ) != 0. else 5 / float(
                min(self.pathogen.rec_rates)
            )
        )
        # set arguments that can be provided but don't need to
        for opt_arg in ['default_susceptibility', 'default_drug_selection_factor']:
            if opt_arg in params:
                setattr(self, '_' + opt_arg, params.pop(opt_arg))
                print '<{0:s}> was specified as an argument.'.format(opt_arg)
        # specify number of decimals to round the log times
        self._log_time_rounding = params.get('log_time_rounding', 2)
        self._ignore_dyn_nodes_in_log = params.get('ignore_dyn_nodes_in_log', False)
        self._time_rounding = params.get('time_rounding', 4)
        self._resolve_hots_pathogen_relations()
        # by default do not consider selection
        self.skip_selection = True
        # run through all the pathogens provided
        for strain_id in xrange(self.pathogen.n):
            # get the rates at which this strain 'mutates' to another strain
            a_selection = self.pathogen.select_rates[strain_id]
            # if this strain changes with a certain rate into another strain
            if any([rate != 0 for rate in a_selection]):
                # we need to consider selection in the simulation
                self.skip_selection = False
        self.treatment = treatment
        # if the treatment argument is provided
        if treatment is None:
            # we can ignore treatment (can always be changed later)
            self.skip_treatment = True
        else:
            # we need to take into account treatment
            self.skip_treatment = False
            self._resolve_treatment_pathogen_relations()
        # initialize the status (-1 means susceptible): everyone is susceptible
        self.current_view = [-1 for _ in xrange(self.contact_structure.n)]
        # indicate whether host was infected through infection or mutation. 0: infection, 1: mutation
        self.current_infection_type = [-1 for _ in xrange(self.contact_structure.n)]  # -1: not infected
        # dict that can be used to pass on values between phases. This is more for future use and not important for now.
        self._phase_passon = {}
        # initialize the status of whether or not a host is under treatment. (-1: not infected
        self.current_treatment = [-1 for _ in xrange(self.contact_structure.n)]
        # this will be a list of booleans (index: strain_id, value: treating yes/no)
        self.treating = []
        self.selecting = []
        # initialize the priorityqueue which will hold all the events (infection, recovering, mutation, ...)
        self.queue = PriorityQueue()
        self._inf_file_o = None
        # This will be the stream to an output file to write incremental steps
        # into.

    @staticmethod
    def _cut_times(recover_time, start_times, stop_times, inf_times, nn):
        """
        Filter the nn list and inf_times for realistic events (within an active connection (dynamic)
            or within the the duration the node is infected (both)
        :param recover_time:
        :param start_times:
        :param stop_times:
        :param inf_times:
        :param nn:
        :return:
        """
        inf_times += start_times  # set the infection times to 'from now'
        # set the end of possible transmission to either the interaction end or the recover time
        stop_times = where(
            stop_times < recover_time,
            stop_times,
            recover_time
        )
        # keep only the neighbours in nn which have an infection time smaller than the nodes recover_time
        return (
            nn[inf_times < stop_times],
            inf_times[inf_times < stop_times],
        )

    def _create_neighbour_events(self, inf_event, nn, inf_times, node_id, token_id):
        """
        This method will be used in any of the (further below defined) event_handler_... methods if we are dealing
            with a static network
        :param inf_event:
        :param nn:
        :param inf_times:
        :param node_id:
        :param token_id:
        :return:
        """
        if not inf_event:  # if this event is a mutation from one Strain to another, update the count
            # of the old Strain
            self.current_infection_type[node_id] = 1  # set the type of infection to mutation
            self.simulation_log['mutations'].append(
                # to do: if we have a dynamic network the degree is replaced with the node_id, ideall we should report
                # the activity of this node
                (self.t, token_id, len(self.contact_structure.nn[node_id])) if self.contact_structure.is_static else
                (self.t, token_id, node_id)
            )
            self._count_per_strains[self.current_view[node_id]] -= 1
        else:
            self.current_infection_type[node_id] = 0  # set the type of infection to 'through selection'
        self.current_view[node_id] = token_id  # set the current status of the node to the Strain he is
        # infected with
        for x in xrange(inf_times.size):  # put all the infection events of neighbours into the queue
            self.queue.put_nowait(Event(self.t + inf_times[x], nn[x], token_id, True, node_id))
        self._count_per_strains[token_id] += 1  # update the counter of the Strain that is infecting.

    def _get_neighbours_static(self, node_id, token_id):
        """
        This method will be used in any of the (further below defined)
        event_handler_... methods if we are dealing with a static network

        Parameter:
        ----------
        :param node_id: the index (id) of the node that is being infected.
        :type node_id: int
        :param token_id: the identifier for the token the node gets.
        :type token_id: int

        :return: tuple with a list of nearest neighbours, the recover time,
            the infection times for each neighbour, the starting time of the
            interactions (here this is 0 as we have a static network) and the
            stop times of each interaction (here this is the recover time
            of the infected node as we are in the static case)

        """
        recover_time = self.pathogen.rec_dists[token_id].get_val()
        nn = []
        inf_times = []
        t_inf = 0
        for n_neigh in self.contact_structure.nn[node_id]:
            t_inf = self.pathogen.trans_dists[token_id].get_val()
            while t_inf <= recover_time:
                nn.append(n_neigh)
                inf_times.append(t_inf)
                t_inf += self.pathogen.trans_dists[token_id].get_val()
        # the last 2 returned arguments are the start and stop condition for
        # the neighbours, which simply boil down to self.t and the recover_time
        # in the static case
        return array(nn), recover_time, array(inf_times), 0., recover_time

    def _get_neighbours_static_single(self, node_id, token_id):
        """
        This method will be used in any of the (further below defined)
        event_handler_... methods if we are dealing with a static network
        and if we constrain a contact between a carrying and a susceptible 
        individual to only ever transmit once.
        Note: You should have a particular reason why to use this function,
            if you don't use the _get_neighbours_static function. You can do so
            by simply not specifying the single_transmission parameter when 
            initializing a scenario.

        Parameter:
        ----------
        :param node_id: the index (id) of the node that is being infected.
        :type node_id: int
        :param token_id: the identifier for the token the node gets.
        :type token_id: int

        :return: tuple with a list of nearest neighbours, the recover time,
            the infection times for each neighbour, the starting time of the
            interactions (here this is 0 as we have a static network) and the
            stop times of each interaction (here this is the recover time
            of the infected node as we are in the static case)

        """
        nn = copy(self.contact_structure.nn[node_id])
        recover_time = self.pathogen.rec_dists[token_id].get_val()
<<<<<<< endemic/nw_spread/Spreading.py
        inf_times = self.pathogen.trans_dists[token_id].v_get(nn) if nn.size \
            else array([])
=======
        inf_times = self.pathogen.trans_dists[token_id].v_get(nn) if nn.size else array([])
        # nn = nn[inf_times < recover_time]
        #inf_times = inf_times[inf_times < recover_time]
        # the last 2 returned arguments are the start and stop condition for the neighbours, which simply boil down to
        # self.t and the recover_time in the static case
>>>>>>> endemic/nw_spread/Spreading.py
        return nn, recover_time, inf_times, 0., recover_time

    def _get_neighbours_dynamic(self, node_id, token_id):
        """
        This method will be used in any of the (further below defined)
        event_handler_... methods if we are dealing with a dynamic network
        :param node_id:
        :param token_id:
        :return:
        """
        recover_time = self.pathogen.rec_dists[token_id].get_val()
        # returns all relevant connections: node_ids, start_times, stop_times
        diff_nn, _start_times, _stop_times = self.contact_structure.get_events(
                node_id,
                self.t,
                recover_time
                )
        # cut the start_times with the current time:
        _start_times = where(
            _start_times >= self.t,
            _start_times - self.t,
            0.0
        )
        _stop_times -= self.t
        nn = []
        inf_times = []
        start_times = []
        stop_times = []
        # TODO: This is list op on np array, not ideal!
        for i in xrange(len(diff_nn)):
            stop_time = _stop_times[i]
            start_time = _start_times[i]
            t_inf = self.pathogen.trans_dists[token_id].get_val()
            while t_inf + start_time <= stop_time:
                nn.append(diff_nn[i])
                inf_times.append(t_inf)
                start_times.append(start_time)
                stop_times.append(stop_time)
                t_inf += self.pathogen.trans_dists[token_id].get_val()
        return array(
                nn
                ), recover_time, array(
                        inf_times
                        ), array(
                                start_times
                                ), array(
                                        stop_times
                                        )

    def _get_neighbours_dynamic_single(self, node_id, token_id):
        """
        This method will be used in any of the (further below defined)
        event_handler_... methods if we are dealing with a dynamic network
        :param node_id:
        :param token_id:
        :return:
        """
        recover_time = self.pathogen.rec_dists[token_id].get_val()
        nn, start_times, stop_times = self.contact_structure.get_events(
                node_id, self.t, recover_time
                )
        # this should return id_list, start_array, stop_array
        inf_times = self.pathogen.trans_dists[token_id].v_get(nn) if nn.size \
            else array([])
        # cut the start_times with the current time:
        start_times = where(
            start_times >= self.t,
            start_times - self.t,
            0.0
        )
        stop_times -= self.t
        return nn, recover_time, inf_times, start_times, stop_times

    # to do: attributes are defined in __init__ so any change in their definitions has to be passed on to here
    # ideally reset should be called in __init__ to avoid any code duplication.
    # ToDo: reset with a new temporal graph is not implemented.
    def reset(self, graph=None):
        """
        This method allow to reset the entire scenario, i.e. the scenario is set back to the time t=start/0 before any
        initial infections took place.
        It can be called when doing multiple simulations with the same parameters.
        If the parameter graph is provided, then the topology of the contact network is updated. (see the method
            update_topology in the ContactNetwork class for details).
        :param graph: Graph object from the nw_construct package
        :return:
        """
        self.log = {}
        self.simulation_log = {
            'phases': {},
            'mutations': [],
            'adjacency': {},
            'modifications': {},
            'param_alternation': {}
        }
        self.outcome = {}
<<<<<<< endemic/nw_spread/Spreading.py
        self.t = self.contact_structure.t_start if not \
                self.contact_structure.is_static else 0
=======
        self.t = self.contact_structure.t_start if not self.contact_structure.is_static else 0
>>>>>>> endemic/nw_spread/Spreading.py
        self._count_per_strains = array([0 for _ in xrange(self.pathogen.n)])
        self._counts_over_time = zeros((1, self.pathogen.n))
        self.current_view = [-1 for _ in xrange(self.contact_structure.n)]  # Indicates the current status of the hosts
        self.current_infection_type = [-1 for _ in xrange(self.contact_structure.n)]
        self.current_treatment = [-1 for _ in xrange(self.contact_structure.n)]
        self._phase_passon = {}  # dict that can be used to pass on values between phases.
        while True:  # empty the queue
            try:
                self.queue.get_nowait()
            except Empty:
                break
        if graph:
            self.contact_structure.update_topology(graph)
        return None

    def _resolve_hots_pathogen_relations(self):
        """
        It might be that some hosts are not or less susceptible to some pathogen strains. This method properly
            established the relation between each host and the present pathogen strains. See ContactNetwork for more
            details.
        :return:
        """
        # run through all the hosts
        for host_id in xrange(len(self.contact_structure._susceptible)):
            # get the susceptibility status for this host
            a_suscept = self.contact_structure._susceptible[host_id]
            # get the default value either from the ContactNetwork object
            if 'Default' in a_suscept:
                default = a_suscept.pop('Default')
            # or from self
            else:
                default = self._default_susceptibility
            # initialize all susceptibilities as the default
            self.contact_structure.susceptible[host_id] = [default for _ in xrange(self.pathogen.n)]
            # if other values are provided (e.g. wild_type: 0.5) convert the pathogen strain name to its id and
            # set the susceptibility for this strain
            for strain_name in a_suscept:
                self.contact_structure.susceptible[host_id][self.pathogen.ids[strain_name]] = a_suscept[strain_name]
        return 0

    def _resolve_treatment_pathogen_relations(self):
        """
        This method defines all treatment related parameters. All relevant parameters are
        stored in self.therapy_* arguments which are in the form of lists. The index of the list
        is equivalent to the therapy id which can be mapped to the therapy name with the
        scenario.treatment.names dict.
        A mapping between pathogen strain and corresponding drug is given by self.therapy_strain_map

        :return:
        """
        nbr_of_therapies = self.treatment.n
        # initialize the probabilities and therapy delays (node specific)
        self.therapy_probas = [[]] * nbr_of_therapies
        self.therapy_delays = [[]] * nbr_of_therapies
        # initialize the list containing for each therapy a factor that will multiply the transmission rate
        self.therapy_trans_facts = [{}] * nbr_of_therapies
        # same thing for the recovery rate
        self.therapy_recover_facts = [{}] * nbr_of_therapies
        # and the same thing for the selection rate (selection rate: rate at which a strain can mutate to another one)
        self.therapy_select_facts = [{}] * nbr_of_therapies
        # the two dicts below will contain a mapping between therapy and strain id and vice versa
        self.therapy_strain_id_map = {}
        self.strain_therapy_id_map = {}
        for a_therapy in self.treatment.therapies:
            # get the therapy id
            its_id = self.treatment.ids[a_therapy.name]
            # initialize the mapping between therapy and strains
            self.therapy_strain_id_map[its_id] = []
            # set the treatment probabilities for each node
            if type(a_therapy.treatment_proba) is float:
                self.therapy_probas[its_id] = [
                    a_therapy.treatment_proba for _ in range(self.contact_structure.n)
                ]
            # to do: implement other cases (not uniform)
            else:
                raise ValueError('Needs to be implemented')
            # self.therapy_probas[its_id] = a_therapy.treatment_proba
            # if the therapy comes with a delay, handle the delay properly
            if type(a_therapy.delay) is float:
                self.therapy_delays[its_id] = [
                    a_therapy.delay for _ in range(self.contact_structure.n)
                ]
            # to do: implement other cases (not uniform)
            else:
                raise ValueError('Needs to be implemented')
            #self.therapy_delays[its_id] = a_therapy.delay
            # self.pathogen.ids is a dict like so: {'wild_type': 0,...} so it links the name to an id
            for strain_name in self.pathogen.ids:
                strain_id = self.pathogen.ids[strain_name]
                # the try except statement is to test whether the strain_name is present in a_therapy.drug.trans-
                #   mission_factor
                try:
                    self.therapy_trans_facts[its_id][
                        strain_id
                    ] = a_therapy.drug.transmission_factor[strain_name]
                    self.therapy_strain_id_map[its_id].append(strain_id)
                # if the strain_name is not present use the default value
                except KeyError:
                    self.therapy_trans_facts[its_id][
                        strain_id
                    ] = a_therapy.drug.transmission_factor['Default']
                # same story for the recover factor
                try:
                    self.therapy_recover_facts[its_id][
                        strain_id
                    ] = a_therapy.drug.recover_factor[strain_name]
                    self.therapy_strain_id_map[its_id].append(strain_id)
                # als here, use default if strain_name is not there
                except KeyError:
                    self.therapy_recover_facts[its_id][
                        strain_id
                    ] = a_therapy.drug.recover_factor['Default']
                # and again same for the selection factor
                try:
                    self.therapy_select_facts[its_id][
                        strain_id
                    ] = a_therapy.drug.selection_factor[strain_name]
                    self.therapy_strain_id_map[its_id].append(strain_id)
                except KeyError:
                    self.therapy_select_facts[its_id][
                        strain_id
                    ] = a_therapy.drug.selection_factor['Default']
                if self.therapy_select_facts[its_id][strain_id] != 1:
                    self.skip_selection = False
                # it might be we added a strain_id several times to self.therapy_strain_id_map[its_id] if so remove
                #   the duplicates
                self.therapy_strain_id_map[its_id] = list(set(self.therapy_strain_id_map[its_id]))
                # initialize the strain therapy map
                self.strain_therapy_id_map[strain_id] = []
        # run through all the therapies
        for therapy_id in self.therapy_strain_id_map:
            # run through all the strain ids for a given therapy (i.e. all the strains that are treated with it)
            for strain_id in self.therapy_strain_id_map[therapy_id]:
                # try to append the therapy id to the mapping in the other direction, i.e. to get all the applied
                # therapies given a strain id
                try:
                    self.strain_therapy_id_map[strain_id].append(therapy_id)
                # if the therapy id was the first one, create the list.
                except KeyError:
                    self.strain_therapy_id_map[strain_id] = [therapy_id]
        return 0

    # just some custom errors to get more specific output if something goes wrong
    class WrongPathogenError(Exception):
        pass

    class NoneImplementationError(Exception):
        pass

    class InitiateInfectionError(Exception):
        pass

    def _initiate_infection(self, strain, ):
        """
        Parameters:
        -----------
        :param strain: dict, key is the name of a strain, value is a list of
                node id's or the name of another pathogen strain or 'random'.
                If the value is 'random' then one random host is infected.
                If the value gives another pathogen strain, then a randomly
                chosen individual that is currently infected with the indicated
                pathogen strain is chosen and will be infected with the strain
                provided in the key.
                If the value is another dict, it can contain the keys 't_inf',
                'host' and 'nbr_infectins'. 't_inf' specifies the time of the
                infection and 'host' can either be a list of node id's, another
                pathogen name or 'random'. The default value of 'host' is
                'random'. 'nbr_infections' specifies how many host should be 
                infected, default is 1.
                Eg. - strain = {'wild_type':[1,5,10]}: infects hosts 1,5 and 10
                        with the wild type strain.
                    - strain = {'mutant_1': 'wild_type'}: infect a randomly
                        chose host that is currently infected with the
                        wild_type strain with the mutant_1 strain.
                    - strain = {'wild_type': {
                            't_inf': 10, 'host': 'random', 'nbr_infection':2
                            }}: 
                        infect two random hosts at time 10.
        """
        for name in strain:
            # make sure the strain exists
            if name not in self.pathogen.ids.keys():
                raise self.WrongPathogenError(
                        "There is no pathogen strain with the name <%s>."%name
                        )
            # if a dictionary is given, infect random node at time t_inf
            if isinstance(strain[name], dict):
                t_inf = strain[name].get(
                        't_inf', self.contact_structure.t_start
                        )
                hosts = strain[name].get(
                        'host', 'random' # defalut is random
                        )
            else:
                t_inf = self.t
                hosts = strain[name]
            # carry out the infection at each t_inf
            def _expander(_keys, _values):
                if isinstance(_keys, list):
                    if isinstance(_values, list):
                        if isinstance(_values[0],list):
                            return _keys, _values
                        else:
                            return _keys, [_values for _ in _keys]
                    else:
                        return _keys, [[_values] for _ in _keys]
                else:
                    if isinstance(_values, list):
                        return [_keys], [_values]
                    else:
                        return [_keys], [[_values]]
            t_inf, hosts = _expander(t_inf, hosts)

            # now the infections are set for the specified strain
            for i in xrange(len(t_inf)):
                a_t_inf = t_inf[i]
                if not self.contact_structure.is_static:
                    candidate_nodes = \
                            self.contact_structure.get_nodes_by_lifetime(
                                a_t_inf
                                )
                else:
                    candidate_nodes = range(self.contact_structure.n)
                if len(candidate_nodes) < 1:
                    raise self.InitiateInfectionError(
                            """

 No host at time %s to be infected.

                            """%a_t_inf
                            )
                # now for all the infection times we have a set of present
                # hosts
                for a_host in hosts[i]:
                    if a_host == 'random':
                        the_host = random.choice(
                                candidate_nodes
                                )
                        self.queue.put_nowait(
                            Event(
                                a_t_inf, the_host, 
                                self.pathogen.ids[name], False,
                                )
                            )
                    # if the host is specified by id
                    elif isinstance(a_host, int):
                        if a_host not in candidate_nodes:
                            raise self.InitiateInfectionError(
                                    """
                                       The host with ID %s does not exist at 
                                       the time it should be infected, %s.
                                    """%(a_host, a_t_inf)
                                    )
                        self.queue.put_nowait(
                                Event(
                                    a_t_inf, a_host, self.pathogen.ids[name],
                                    False,
                                    )
                                )
                    # the host must be a specific strain
                    else:
                        if a_t_inf != self.t:
                            raise self.InitiateInfectionError(
                                    """
                                        The targeted infection of a host 
                                        infected with a specific pathogen
                                        is only possible if the infection 
                                        time is the current time of the 
                                        simulation.
                                        Current time: %s
                                        Time of infection: %s
                                    """%(self.t, a_t_inf)
                            )
                        # get the ids of of name and hosts
                        target_id = self.pathogen.ids[a_host]
                        mut_id = self.pathogen.ids[name]
                        # check if the target pathogen is present in the
                        # population, if not raise an error.
                        potentials = [
                            filter(
                                lambda x: x in filter(
                                    lambda x: self.current_view[x] == target_id,
                                    range(self.contact_structure.n)
                                ), candidate_nodes
                            )
                        ]
                        if not potentials:
                            raise self.InitiateInfectionError(
                                    """
                                        There are no host infected with %s at 
                                        the moment.
                                    """%a_host
                                    )
                        # chose at random one specified infected
                        # host -> new_mutated.
                        new_mutated = random.choice(potentials)
                        # now that we have chose our host to mutate, we need to
                        # filter the existing event queue for any events
                        # associated with this host (e.g. this host will not
                        # infect its neighbours with its previous strain
                        # anylonger.
                        # copy the event queue into a list
                        pending_events = []
                        while True:
                            try:
                                pending_events.append(self.queue.get_nowait())
                            except Empty:
                                break
                        # remove all entries where new_mutated infects plus the
                        # one where he recovers an recreate the queue
                        for an_event in pending_events:
                            # filter out infections where new_mutated is the
                            # source
                            if an_event[1][3] != new_mutated:
                                #do not take the recover event for new_mutated
                                if an_event[1][0] != new_mutated:  
                                    self.queue.put_nowait(an_event)
                        # add infection event of new_mutated with hosts
                        # make sure that the infection overwrites new_mutated's
                        # old status (use inf_event = False)
                        self.queue.put_nowait(
                                Event(a_t_inf, new_mutated, mut_id, False, )
                                )
        return 0

    # # this method is not used anymore. Could be removed.
    # def initiate_infection(self, strain, ):
    #     """
    #     Function to set the initial seed for an infection.
    #     
    #     Arguments:
    #         - strains: dict, key the name of a strain, value a list of node id's or
    #             'random'. If the value is 'random' then one random host is infected.
    #             Eg. strain = {'wild_type':[1,5,10]}: infects _hosts 1,5 and 10 with the
    #                 wild type strain.
    #     """
    #     self.t = 0
    #     for name in strain:
    #         if name not in self.pathogen.ids.keys():
    #             raise self.WrongPathogenError("There is no pathogen strain with the name <%s>." % name)
    #         if type(strain[name]) is not str:
    #             for node_id in strain[name]:
    #                 self.current_view[node_id] = self.pathogen.ids[name]
    #         else:
    #             self.current_view[nrand.randint(0, self.contact_structure.n)] = self.pathogen.ids[name]
    #         self._init_queue()
    #     return 0

    # # unused method can be removed (along with self.initiate_infection)
    # def _init_queue(self, ):
    #     """
    #     Initiate the priority queue according to self.current_view
    #     """
    #     for node_id in xrange(self.contact_structure.n):
    #         if self.current_view[node_id] != -1:
    #             self.queue.put_nowait(Event(self.t, node_id, self.current_view[node_id], True, ))
    #             self.current_view[node_id] = -1
    #     return 0

    # here below follow several _handle_event... functions each one of these take an event (node id, token id, inf type,
    # source) as an argumnet (see Event class for further details) and digest it. Based on the event self.current_view
    # will be updated an new events will be created

    # method to handle events if we have both treatment and mutation/selection. Maybe best start with the method
    # self._handle_event_simple as this is for the most trivial case (no treatment, no selection/mutation)
    def _handle_event_combined(self, an_event, get_neighbours):  # with selection & treatment
        """
        This method handles the events in a spreading process with selection (mutation + selection) and treatment.

        :param an_event:
        :return:
        """
        node_id, token_id, inf_event, source = an_event
        if token_id == -1:  #the Event is recovering
            old_strain_id = self.current_view[node_id]
            self.contact_structure.susceptible[node_id][old_strain_id] = self.pathogen.rec_types[old_strain_id]
            # set the node status back to susceptible
            self.current_view[node_id] = -1
            self.current_infection_type[node_id] = -1
            self._count_per_strains[old_strain_id] -= 1
        else:  #the Event is an infection or selection
            if inf_event and self.current_view[node_id] != -1:  #infection of infected host: do nothing
                pass  #NOTE: if super-infections are possible, here they would take place
            else:  #infection of susceptible host or mutation
                if nrand.rand() < self.contact_structure.susceptible[node_id][token_id] or not inf_event:
                    nn, recover_time, inf_times, start_times, stop_times = get_neighbours(node_id, token_id)
                    # This is the part devoted to selection and treatment
                    # ##
                    if self.selecting[token_id]:
                        selected_strain_id, selection_times = self.pathogen.get_selected(token_id)
                        if selection_times[selected_strain_id] < recover_time:
                            recover_time = selection_times[selected_strain_id]
                            new_token, new_inf_event = selected_strain_id, False
                        else:
                            new_token, new_inf_event = -1, True
                    else:
                        new_token, new_inf_event = -1, True

                    if self.treating[token_id]:
                        therapy_ids = self.strain_therapy_id_map[token_id]
                        #to do: gather the various times and chose at random one, not the
                        # smallest as now.
                        #print self.therapy_select_facts
                        #issue: this does not work if we have more than one therapy.
                        for therapy_id in therapy_ids:
                            if nrand.rand() < self.therapy_probas[therapy_id][node_id]:
                                delay = self.therapy_delays[therapy_id][node_id]
                                if recover_time > delay:  #will recover after treatment delay
                                    recover_time = delay + (
                                                               recover_time - delay
                                                           ) * self.therapy_recover_facts[therapy_id][token_id] ** (-1)
                                    if self.selecting[token_id]:
                                        selection_times = [
                                            delay +
                                            (selection_times[x] - delay) *
                                            self.therapy_select_facts[therapy_id][x] ** (-1)
                                            for x in xrange(len(selection_times))
                                        ]  #x is the id of the potential mutant
                                        selected_strain_id = selection_times.index(min(selection_times))
                                        if recover_time > selection_times[selected_strain_id]:
                                            recover_time = selection_times[selected_strain_id]
                                            new_token, new_inf_event = selected_strain_id, False
                                inf_times = where(
                                    start_times + inf_times <= delay,
                                    inf_times,
                                    delay + (inf_times - delay) * self.therapy_trans_facts[therapy_id][token_id] ** (-1)
                                )
                                # ##
                    nn, inf_times = self._cut_times(recover_time, start_times, stop_times, inf_times, nn)
                    self.queue.put_nowait(Event(self.t + recover_time, node_id, new_token, new_inf_event, ))
                    self._create_neighbour_events(inf_event, nn, inf_times, node_id, token_id)
        return 0

    # no selection, no treatment (but might still need to handle selections)
    def _handle_event_simple(self, an_event, get_neighbours):
        """
        This method handles events in a spreading process without treatment nor selection.

        :param an_event:
        :return:
        """
        node_id, token_id, inf_event, source = an_event  # token_id is the id of the pathogen -1 means recovering/dead
        if token_id == -1:  # the Event is recovering
            old_strain_id = self.current_view[node_id]  # what was the old status of that node
            # self.pathogen.rec_types: a list indicating how one recovers after an infection. The index is the pathogen
            # id and the value is either 0,1 meaning ether back to susceptible or resistant.
            self.contact_structure.susceptible[node_id][old_strain_id] = self.pathogen.rec_types[old_strain_id]
            self.current_view[node_id] = -1  # set the node back to the uninfected state
            self.current_infection_type[node_id] = -1  # set the infection type back
            self._count_per_strains[old_strain_id] -= 1  # update the count of number of infected for that strain
        else:  # the Event is an infection
            if inf_event and self.current_view[node_id] != -1:  # infection of infected host: do nothing
                pass  # NOTE: if super-infections are possible, here they would take place
            else:  # infection of susceptible host or mutation
                # this reads: if the node is susceptible (this is all before the 'or not') or it is actually not an
                # infection event (a mutation in this case)
                if nrand.rand() < self.contact_structure.susceptible[node_id][token_id] or not inf_event:
                    nn, recover_time, inf_times, start_times, stop_times = get_neighbours(node_id, token_id)
                    # This is the method without selection nor treatment, so not much to be done here
                    # ##
                    # ##
                    nn, inf_times = self._cut_times(recover_time, start_times, stop_times, inf_times, nn)
                    self.queue.put_nowait(Event(self.t + recover_time, node_id, -1, True, ))  # put the recover event
                    self._create_neighbour_events(inf_event, nn, inf_times, node_id, token_id)
        return 0

<<<<<<< endemic/nw_spread/Spreading.py

    # this is the incremental version of the simple event handler. In this
    # version each change is written to an output file
    def _handle_event_simple_inc(self, an_event, get_neighbours):
        """
        This method handles events in a spreading process without treatment nor
        selection.

        :param an_event:
        :param get_neighbours:
        :param inc_file: opened file to write changes into
        :return:
        """
        node_id, token_id, inf_event, source = an_event
        # token_id is the id of the pathogen -1 means recovering/dead
        if token_id == -1:
            # the Event is recovering
            old_strain_id = self.current_view[node_id] 
            # what was the old status of that node.
            # self.pathogen.rec_types: a list indicating how one recovers
            # after an infection. The index is the pathogen id and the value is
            # either 0,1 meaning ether back to susceptible or resistant.
            self.contact_structure.susceptible[
                    node_id
                    ][old_strain_id] = self.pathogen.rec_types[old_strain_id]
            self.current_view[node_id] = -1 
            # set the node back to the uninfected state
            self.current_infection_type[node_id] = -1 
            # set the infection type back
            self._count_per_strains[old_strain_id] -= 1 
            # update the count of number of infected for that strain
            self._inc_file_o.write(
                    '%s, %s\n' % (
                        self.t, self.contact_structure.all_nodes[node_id]
                        ) 
                    )
        else:
            # the Event is an infection
            if inf_event and self.current_view[node_id] != -1: 
                # infection of infected host: do nothing
                pass 
                # NOTE: if super-infections are possible, here they would take
                # place
            else: 
                # infection of susceptible host or mutation
                if nrand.rand() < self.contact_structure.susceptible[node_id][
                        token_id
                        ] or not inf_event:
                    nn, recover_time, inf_times, start_times, stop_times = \
                            get_neighbours(node_id, token_id)
                    # : if the node is susceptible (this is all before the 'or
                    # not') or it is actually not an infection event (a
                    # mutation in this case).
                    self._inc_file_o.write(
                            '%s, %s, %s\n' % (
                                self.t,
                                self.contact_structure.all_nodes[node_id],
                                self.contact_structure.all_nodes[source] if
                                source is not None else 'seed'
                                )
                            )
                    # This is the method without selection nor treatment, so
                    # not much to be done here
                    nn, inf_times = self._cut_times(
                            recover_time, start_times, stop_times, inf_times, nn
                            )
                    self.queue.put_nowait(
                            Event(self.t + recover_time, node_id, -1, True,)
                            )
                    # put the recover event
                    self._create_neighbour_events(
                            inf_event, nn, inf_times, node_id, token_id
                            )
                    # cerate and add the infection events for the neighbours.
        return 0

    def _handle_event_selection(self, an_event, get_neighbours):  #with only selection
=======
    def _handle_event_selection(self, an_event, get_neighbours):  # with only selection
>>>>>>> endemic/nw_spread/Spreading.py
        """
        This method handles events in a spreading process with selection (mutation + selection) but without treatment.

        :param an_event:
        :return:
        """
        node_id, token_id, inf_event, source = an_event
        if token_id == -1:  #the Event is recovering
            old_strain_id = self.current_view[node_id]
            self.contact_structure.susceptible[node_id][old_strain_id] = self.pathogen.rec_types[old_strain_id]
            self.current_view[node_id] = -1
            self.current_infection_type[node_id] = -1
            self._count_per_strains[old_strain_id] -= 1
        else:  #the event is an infection or selection
            if inf_event and self.current_view[node_id] != -1:  #infection of infected host: do nothing
                pass  # NOTE: if super-infections are possible, here they would take place
            else:  #infection of susceptible host or selection/mutation
                if nrand.rand() < self.contact_structure.susceptible[node_id][token_id] or not inf_event:
                    nn, recover_time, inf_times, start_times, stop_times = get_neighbours(node_id, token_id)
                    # This is the part devoted to selection
                    # ##
                    # determine the strain that is selected for and the time at which the mutation will take place
                    # see Pathogen.get_selected method for more details
                    new_token, new_inf_event = -1, True  # if the mutation arises after recovering, the subsequent
                    # event is simply: recovered
                    # infections of the neighbours is now as without the selection/mutation as we assured that recover
                    # time is either the true recover time or the mutation time.
                    # ##
                    if self.selecting[token_id]:
                        selected_strain_id, selection_times = self.pathogen.get_selected(token_id)
                        if selection_times[
                            selected_strain_id] < recover_time:  #if the mutation is before the recovering
                            recover_time = selection_times[
                                selected_strain_id]  # adjust the time of "recover" from the current infection.
                            new_token, new_inf_event = selected_strain_id, False  #set the token and infection event status
                            # for a subsequent event.
                    nn, inf_times = self._cut_times(recover_time, start_times, stop_times, inf_times, nn)
                    self.queue.put_nowait(Event(self.t + recover_time, node_id, new_token, new_inf_event, ))
                    # when writing new_token and new_inf_event into the queue, it is either just the recover event
                    # or the mutation event.
                    self._create_neighbour_events(inf_event, nn, inf_times, node_id, token_id)
        return 0

    def _handle_event_treatment(self, an_event, get_neighbours):  # with only treatment
        """
        This method handles the events in a spreading process with treatment but without selection.

        :param an_event:
        :return:
        """
        node_id, token_id, inf_event, source = an_event
        if token_id == -1:  # the Event is recovering
            old_strain_id = self.current_view[node_id]
            # issue: import those into the namespace of spreading: eg. self._contact_network__susceptible...
            self.contact_structure.susceptible[node_id][old_strain_id] = self.pathogen.rec_types[old_strain_id]
            self.current_view[node_id] = -1
            self.current_infection_type[node_id] = -1
            self._count_per_strains[old_strain_id] -= 1
        else:  # the Event is an infection or selection
            if inf_event and self.current_view[node_id] != -1:  # infection of infected host: do nothing
                pass  # NOTE: if super-infections are possible, here they would take place
            else:  # infection of susceptible host or mutation/selection
                if nrand.rand() < self.contact_structure.susceptible[node_id][token_id] or not inf_event:
                    nn, recover_time, inf_times, start_times, stop_times = get_neighbours(node_id, token_id)
                    # This is the part devoted to treatment
                    # ##
                    therapy_ids = self.strain_therapy_id_map[token_id]  # knowing the id of the pathogen strain, we
                    # get the id of the therapy (or therapies) that applies to this strain.
                    # issue: at the moment this approach takes only into account one therapy per pathogen strain
                    for therapy_id in therapy_ids:
                        # if this pathogen is treated at the moment and in this particular case we treat:
                        # To Do: see combined case: you can move self.treating condition outside the therapy_id loop
                        if self.treating[token_id] and nrand.rand() < self.therapy_probas[therapy_id][node_id]:
                            delay = self.therapy_delays[therapy_id][node_id]  # get the delay until treatment
                            if recover_time > delay:  # if node has not recovered during the delay
                                # # define a new recover time, that is the delay plus a new recover time under treatment
                                # print 'rec_1', recover_time
                                recover_time = delay + (recover_time - delay) * \
                                                       self.therapy_recover_facts[therapy_id][token_id] ** (-1)
                                # print 'rec_2', recover_time
                                # # define potential new infection times
                                # new_inf_times = delay + \
                                # self.pathogen.trans_dists[token_id].v_get(nn) if nn.size else array([])
                                # if the infection time is bigger than the delay, take the modified infection time.
                                # print 'therapy factor', self.therapy_trans_facts[therapy_id][token_id]
                                # print 'inf_1', inf_times
                                inf_times = where(
                                    start_times + inf_times <= delay,  # from now until infection is smaller
                                    inf_times,
                                    # delay + new_inf_times
                                    delay + (inf_times - delay) * self.therapy_trans_facts[therapy_id][token_id] ** (-1)
                                )
                                # get the time from now until start and stop
                    # ##
                    nn, inf_times = self._cut_times(recover_time, start_times, stop_times, inf_times, nn)
                    self.queue.put_nowait(Event(self.t + recover_time, node_id, -1, True, ))
                    self._create_neighbour_events(inf_event, nn, inf_times, node_id, token_id)
        return 0

    # this method makes things moving.
    def run(self, phases):
        """
        This method will work through a list of phases. Each phase is a dictionary and each entry in this dictionary can
        be seen as a task. You might just specify a single task for a phase or several but in order to keep things
        tractable, use as few tasks as possible for a phase.
        Below are all possible tasks listed.

        Possible tasks:

            'new_infection': {strain_name: 'random'/IDs/other_strain}
                This will introduce a strain into the population. the strain_name specifies which strain should be
                introduced, 'random'/IDs/other_strain will specifies who should be infected.
                (see self._initiate_infection method for details)
            'parameter_alternation': {strain_name: {some_rate: value}}}
                This will change parameters of the pathogen.
                (see self._parameter_alternation method for detail
                Note: if Default values were set for the mutation selection rate they need to be set again.
            't_start': float
            't_stop': float
                If t_stop is provided the method will try to take out events from the event queue until it finds an
                event with t <= t_stop (or some other halt condition is met - like the assert_survival task below)
                Note: it might be that a phase does not run until t_stop but ends earlyer (e.g. if 'builing_up' task
                if provided - see below). For such cases we cannot tell at what time the phase ends and thus at what
                time the next phase should start. For this there is the self._pase_passon attribute which will be used

            'dt': float
                Can be used to specify the time interval for which the current status should be reported. This only
                really is usefull if 'explicit': True (see below).
            'with_treatment': boolean
            'with_selection': boolean
                If provided will force the scenario to use the specific event handler. If not provided both default to
                False.
            'treating': [strain_name, ...]
                If only some strains should be treated you can specify them with the 'treating' task. If 'treating' is
                not provided but 'with_treatment': True, all strains are treated (if the scenario has drug for them).
                So, best always provide the 'treating' task.
            'assert_survival': [strain_name, ...]
                If this task is provided, the simulation stops if one of the specified strain dies out.
            'halt_condition': [strain_name, ...]
                If provided the halt_condition task will check if the specified strains are in QSS and if so, stop the
                simulation.
            'building_up': {strain_name: value, ref} the value indicates the amount of refs that should be infected with
                the indicated strain. If during the phase this values is attained (or exceeded) the phase is stopped and
                the scenario time (self.t) is returned.
                The building up can either be relative or absolute.
                    For relative building up, a typical configuration could look as follows:
                        'building_up': {'relative': {mutant_1: (0.9, total_infected)}}
                        This reads: the building up stops once the mutant_1 strain reached 90% of the total number of
                        infected. Alternatively 'total_infected' can be replaced by 'total' which would mean the total
                        number of hosts or by any other name of present strains, e.g. the wild_type.
            'explicit': boolean.
                This task forces the scenario to provide detailed status reports over the entire simulation.
                if True, then on every self.dt the current status is written into self.log (slowdown!).
            'incremental': string.
                This task only works if explicit==True. It will write every
                event to the specified output file. Note that the file is
                opened in append mode so to be able to write several phases
                into the same file. Be sure to change the file name for each
                run as otherwise several runs might be written into the same
                file. 
            'shuffle': dict. with 'source', 'target', 'substitute' and 'mode':
                This task will shuffle in the specified way the infection status of hosts.
                'target' must be a list of pathogen names and/or 'susc' for susceptible indicating the group of hosts
                    withing which a shuffle should take place.
                'source' a list of pathogen names and/or 'susc' for susceptible indicating the group that should be
                    distributed within the target group.
                'substitute' a list of pathogen names and/or 'susc' if this list is not empty, all hosts of the source
                    group will get assigned randomly a state from the substitute list.
                'mode': optional argument (default value is 'keep') determines how the recover times should be handled.
                    possible values are 'keep', 'reset', 'source', 'both'.
                        'keep' will not reset the recover times and 'reset' will reset all the recover times. 'source'
                        only resets the recover times of all the source nodes and 'both' will reset the recover times
                        for both 'source' and 'target'
                example: 'shuffle': {'target' = [wild_type, mutant_1], 'source': [wild_type, mutant_1],
                    'substitute': []}
                    This will redistribute the status of all the hosts infected with either wild_type or mutant_1 among
                    themselves.
            'state_change': dict with either 'node' or 'relative' and an optional 'mode'. If 'node', for each node id
                give a new token. new tokens are either pathogen names or 'susc' for susceptible (later 'r' for
                    resistant)
                If 'relative' a dict needs to be provided, specifying what token to replace (key) with how much of what
                    e.g. 'wild_type': (0.5, 'mutant_1') reads replace 0.5 of the wild_types with a mutant_1 token.
                If 'mode' is present the event queue is reset with the given mode. If mode is 'keep' then the newly
                introduced tokens will keep the recover times of their previous token, with exception of those that were
                initially susceptible. if mode is 'reset' all the recover times are reset.
                Infection times for neighbouring nodes are reset no matter what.
                example: - 'state_change': {1: 'mutant_1', 2: 'wild_type', ...} will simply change the token of node 1
                        to 'mutant_1' and of node 2 to 'wild_type'. Note that the pathogenic parameters are not reset,
                        i.e.  node 1 will keep its previous recover time and continue to infect with its previous token.
                         - 'state_change': {'mode': 'reset', 1: 'mutant_1'} will change the token for node 1 to mutant_1
                        and draw a new recovery time for node_1 as well as new infection events for its neighbours.
            'reset_recovery': boolean. If true, all the recover times will be reset (drawn from the recover rate dists)
            'add_treatment': adds a new treatment to the existing scenario. Note that the treatment status is not
                changed, i.e. just adding a new treatment will not lead to an actual treatment. You'll need to provide
                the 'with_treatment': True task to run the simulation with a treatment.
        """
        self.simulation_log['scenario'] = phases
        self.simulation_log['adjacency'][self.t] = self.contact_structure.nn
        # issue: should the adjacency be written in each phase?
        self._update_phase_in_sim_log()
        for phase in deepcopy(phases):
            with_run = True  # will call self._run for this phase
            # check if previous phase imposes parameters:
            if self._phase_passon:
                to_pass = self._phase_passon.keys()
                for param in to_pass:
                    phase[param] = self._phase_passon.pop(param)
            if 'parameter_alternation' in phase:
                # note: only transmission_rate and recover_rate changes work.
                alternations = phase.pop('parameter_alternation')
                self.simulation_log['param_alternation'] = alternations
                self._parameter_alternation(alternations)
                with_run = False
            if 'new_infection' in phase:
                infection = phase.pop('new_infection')
                self._initiate_infection(infection)
                try:
                    self.simulation_log['modifications']['new_infection'] = (self.t, infection)
                except KeyError:
                    self.simulation_log['modifications'] = {'new_infection': (self.t, infection)}
                with_run = False
            if 'introducing' in phase:
                # to do: introduce with a certain rate
                introduction_scenario = phase.pop('introducing')
                with_run = False
                # which: which strain is introduced
                # rate: at which rate it is introduced
                # target: set of potential target hosts
            if 'reset_recovery' in phase:
                if phase['reset_recovery']:
                    self._sync_event_queue(mode='reset')
                with_run = False
            if 'state_change' in phase:
                state_change = phase.pop('state_change')
                try:
                    self.simulation_log['modifications']['state_change'] = (self.t, state_change)
                except KeyError:
                    self.simulation_log['modifications'] = {'state_change': (self.t, state_change)}
                # define the mode that is later passed to the sync_event_queue fct
                if 'mode' in state_change:
                    to_mode = state_change.pop('mode')
                else:
                    to_mode = 'changed'
                for a_key in state_change:
                    if type(a_key) is int:
                        self.current_view[a_key] = self.pathogen.ids[state_change[a_key]]
                        try:
                            nodes_to_sync[a_key] = 1
                        except NameError:
                            nodes_to_sync = [0] * self.contact_structure.n
                            nodes_to_sync[a_key] = 1
                        if to_mode == 'changed':
                            to_mode = nodes_to_sync
                    elif a_key == 'relative':
                        rel_replacement = state_change['relative']
                        # 'wt': (0.5, 'mt')
                        for to_replace in rel_replacement:
                            old_token = self.pathogen.ids[to_replace] if to_replace in self.pathogen.ids else -1
                            fraction_to_replace, new_token_name = rel_replacement[to_replace]
                            # the new status is either a pathogen or susceptible
                            new_token = self.pathogen.ids[new_token_name] if new_token_name in self.pathogen.ids else -1
                            nodes_old_token = filter(
                                lambda i: self.current_view[i] == old_token, range(len(self.current_view))
                            )
                            if not len(nodes_old_token):
                                raise ValueError('No nodes are present with the old token')
                            # choose a fraction_to_replace nodes among the nodes_old_token
                            length = len(nodes_old_token)
                            # always round down to the next int
                            nbr_to_replace = int(fraction_to_replace * length)
                            # replace the nodes token
                            for i in xrange(nbr_to_replace):
                                self.current_view[nodes_old_token[i]] = new_token
                                try:
                                    nodes_to_sync[nodes_old_token[i]] = 1
                                except NameError:
                                    nodes_to_sync = [0] * self.contact_structure.n
                                    nodes_to_sync[nodes_old_token[i]] = 1
                        if to_mode == 'changed':
                            to_mode = nodes_to_sync
                    else:
                        raise ValueError(
                            'The state_change phase is not understood, please either use "relative" or node ids.'
                        )
                # update the strain counts
                for strain_id in self.pathogen.ids.values():
                    self._count_per_strains[strain_id] = self.current_view.count(strain_id)
                # reset the event queue
                if to_mode == 'changed':  # if the mode has not been adjusted yet (this should not happen)
                    to_mode = 'keep'
                self._sync_event_queue(mode=to_mode)
                with_run = False
            if 'shuffle' in phase:
                shuffle = phase.pop('shuffle')
                if 'mode' in shuffle:
                    mode = shuffle.pop('mode')
                else:
                    mode = 'keep'
                # write it to the simulation log
                # print self.current_view
                #print self.current_view.count(0), self.current_view.count(1)
                #print self._count_per_strains
                try:
                    self.simulation_log['modifications']['shuffle'] = (self.t, shuffle)
                except KeyError:
                    self.simulation_log['modifications'] = {'shuffle': (self.t, shuffle)}
                # get the id for target and source tokens
                target_ids = [self.pathogen.ids[name] if name != 'susc' else -1 for name in shuffle['target']]
                source_ids = [self.pathogen.ids[name] if name != 'susc' else -1 for name in shuffle['source']]
                # get the substitute token (if any)
                if 'substitute' in shuffle and shuffle['substitute']:
                    substitute_ids = [
                        self.pathogen.ids[name] if name != 'susc' else -1 for name in shuffle['substitute']
                    ]
                else:
                    substitute_ids = []
                # lists of all the nodes with target resp source tokens
                targets = []
                sources = []
                # list of all the tokens to redistribute
                sources_tokens = []
                # build up the new current view
                new_current_view = []
                current_view = self.current_view
                for node in xrange(len(current_view)):
                    new_current_view.append(current_view[node])  # this might change later
                    if current_view[node] in target_ids:
                        targets.append(node)
                    if current_view[node] in source_ids:
                        sources.append(node)
                        # add the node's token to the token to redistribute list
                        sources_tokens.append(current_view[node])
                        # if we have substitutes for the source token, choose a random substitute
                        if substitute_ids:
                            new_current_view[-1] = _get_rand_el(substitute_ids)
                # as long as not all token have been redistributed
                while len(sources_tokens):
                    # pick a node from the targets
                    a_node = _pick_rand_el(targets)
                    # get this node a new token
                    new_token = _pick_rand_el(sources_tokens)
                    # give the node the new token
                    new_current_view[a_node] = new_token
                # the token are redistributed now
                # update the current view with the new view
                self.current_view = new_current_view
                # recreate the event queue (adding new infection events)
                # to do: the recreation of the queue works only without treatment and selection here
                if mode in ['keep', 'reset']:
                    self._sync_event_queue(mode=mode)
                else:
                    raise self.NoneImplementationError(
                        'shuffle phase with a mode different from "keep" or "reset" has not yet been implemented.'
                    )
                #print self.current_view
                # update the strain counts
                for strain_id in self.pathogen.ids.values():
                    self._count_per_strains[strain_id] = self.current_view.count(strain_id)
                with_run = False
            if 'add_treatment' in phase:
                treatment_to_add = phase.pop('add_treatment')
                if self.treatment is not None:
                    self.treatment.add_new_treatment(treatment_to_add)
                else:
                    self.treatment = treatment_to_add
                    self.skip_treatment = False
                self._resolve_treatment_pathogen_relations()
            # most of the specified tasks are taken care of by now. However, if t_stop is provided we still need to
            # carry out the actual simulation. This is done with the self._run method
            if with_run:
                t_start = phase.pop('t_start', None)
                t_stop = phase.pop('t_stop', None)
                # counts_length = len(self._counts_over_time)
                # if counts_length < int(t_stop):
                # self._counts_over_time.extend([None for _ in xrange(int(t_stop) - len(self._counts_over_time))])
                #    for _ in xrange(max(0, int(t_stop) - counts_length) + 10):  #+10 is just a margin
                #        self._counts_over_time = n_append(self._counts_over_time, [[0 for _ in xrange(self.pathogen.n)]], 0)

                # call the self._run method and pass the remaining tasks on to it
                self._run(
                    t_start=t_start,
                    t_stop=t_stop,
                    **phase
                )
                if self._inf_file_o:
                    self._inf_file_o.close()
                    # close the file stream to the incremental output file if
                    # it exists.
                    self._inf_file_o = None
            self._update_phase_in_sim_log()
            self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
            try:
                self.outcome[self.t].append(self.get_outcome)
            except (KeyError, AttributeError):
                self.outcome[self.t] = [self.get_outcome]
                # if 'assert_survival' in phase:
                # to_survive = phase.pop('assert_survival')
                #    break_here = False
                #    for strain_name in to_survive:
                #        if strain_name not in self.outcome['logs'][self.t]['abundance'].keys():
                #            break_here = True
                #            print '%s got extinct.'%(strain_name)
                #    if break_here:
                #        break
        return 0

    def _run(self,
             t_stop,
             with_selection=False,
             with_treatment=False,
             halt_condition=False,
             assert_survival=False,
             t_start=None,
             **params
    ):
        """
        Run the scenario from t_start to t_stop with the given conditions

        :param assert_survival: If list of strain names, it is checked for survival of those strains.
        :param t_start: Gives the starting time
        :param t_stop: Gives the end time.
        :param with_selection: If True, the selection/mutation is considered.
        :param with_treatment: If True, the specified treatment is applied.
            NOTE: If True, consider to put the <treating> argument to specify which strains are treated.
        :param halt_condition: If list of strain names, it is checked for quasi_stability for those strains
        :param params: Collection of optional arguments (remaining tasks):
            - treating: Dict indicating which strain is treated, eg. treating={'wild_type': True, 'Default': False}
                If 'Default' is given, this value will be applied to all strains missing in the dictionary.
            - dt: The time interval after which to update self.log
        """
        # self.t = round(t_start, 4)
        #specify the new run in self.simulation_log
        #self._update_phase_in_sim_log(
        #    selection=with_selection,
        #    therapy=with_treatment,
        #    t_stop=t_stop,
        #    assert_survival=assert_survival,
        #    params=params
        #)
        # determine if the network is static or dynamic and set the appropriate methods.
        if self.contact_structure.is_static:
            if self.single_transmission:
                get_neighbours = self._get_neighbours_static_single
            else:
                get_neighbours = self._get_neighbours_static
        else:
            if self.single_transmission:
                get_neighbours = self._get_neighbours_dynamic_single
            else:
                get_neighbours = self._get_neighbours_dynamic
        # define the event_handler as the simple method for now. This will be adapted if needed in the next lines
        if 'incremental' in params:
            event_handler = self._handle_event_simple_inc
            self._inc_file_o = open(params['incremental'], 'a')
        else:
            event_handler = self._handle_event_simple
        self.treating = []
        self.selecting = []
        # if no selection parameters where provided when initializing the scenario no selection will be attempted no
        # matter what is specified in the phase we are currently in.
        #if self.skip_selection:
        #    with_selection = False
        # same goes for treatment. If no treatment was specified when initializing (and the task 'add_treatment' was
        # never provided so far) we skip the treatment, no mather what is specified in the current phase.
        if self.skip_treatment:
            with_treatment = False
        if with_treatment:
            # if we have treatment, clarify which strain to treat
            if 'treating' in params:
                treat_dict = params.pop('treating')
                def_val = True
                if 'Default' in treat_dict:
                    def_val = treat_dict.pop('Default')
                self.treating = [def_val for _ in xrange(self.pathogen.n)]
                for strain_name in treat_dict:
                    self.treating[self.pathogen.ids[strain_name]] = treat_dict[strain_name]
            else:
                self.treating = [True for _ in xrange(self.pathogen.n)]  #if treating is missing all strains are treated
            # if we have treatment and selection, we need to use the combined event handler
            if with_selection:
                if 'selecting' in params:
                    selecting_dict = params.pop('selecting')
                    def_val = True
                    if 'Default' in selecting_dict:
                        def_val = selecting_dict.pop('Default')
                    self.selecting = [def_val for _ in xrange(self.pathogen.n)]
                    for strain_name in selecting_dict:
                        self.selecting[self.pathogen.ids[strain_name]] = selecting_dict[strain_name]
                else:
                    self.selecting = [True for _ in xrange(self.pathogen.n)]
                event_handler = self._handle_event_combined
            # if it is only treatment, the treatment event handler is the one to use
            else:
                event_handler = self._handle_event_treatment
        elif with_selection:
            if 'selecting' in params:
                selecting_dict = params.pop('selecting')
                def_val = True
                if 'Default' in selecting_dict:
                    def_val = selecting_dict.pop('Default')
                self.selecting = [def_val for _ in xrange(self.pathogen.n)]
                for strain_name in selecting_dict:
                    self.selecting[self.pathogen.ids[strain_name]] = selecting_dict[strain_name]
            else:
                self.selecting = [True for _ in xrange(self.pathogen.n)]
            # at this point we know that only selection is on, so use the selection event handler.
            event_handler = self._handle_event_selection
        # check if the time interval for reporting is specified, if not use default one.
        if 'dt' in params:
            dt = params.pop('dt')
        else:
            dt = self._dt
        # get the next time at which the status of the hosts should be reported or halt condition to be checked
        t_next_bin = self.t + dt
        # check if we have some halt conditions that might lead to a termination of the simulation before t_stop.
        with_halt_condition = False
        if halt_condition or assert_survival:
            with_halt_condition = True
            focus_strain_ids = None
            surviving_strain_ids = None
        if halt_condition:
            focus_strain_ids = array([self.pathogen.ids[strain_name] for strain_name in halt_condition])
        if assert_survival:  # stop as soon as one of the specified stains goes extinct
            surviving_strain_ids = array([self.pathogen.ids[strain_name] for strain_name in assert_survival])
            # with_logging = params.get('explicit', False)
            done = False

            # TO DO: start for new structure.The wile loop can be put after the running conditions and
            # each condition defines its proper stepper function.
            # TO DO: problem of combining conditions remains.

            def stepper(self):
                # get the next event
                (time, n_event) = self.queue.get_nowait()
                # update the time of the scenario
                self.t = round(time, self._time_rounding)
                #self._counts_over_time[int(self.t)] = self._count_per_strains
                # pass the event to the event handler
                event_handler(n_event, get_neighbours)
                # the new time is after the checking time
                if self.t >= t_next_bin:
                    # check for the condition
                    for strain_id in surviving_strain_ids:
                        if not self.current_view.count(strain_id):
                            return 1
                return 0

            while self.t < t_stop:
                try:
                    if stepper(self):
                        break
                except Empty:
                    self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                    break
            """
            while self.t < t_stop and not done:
                try:
                    # get the next event
                    (time, n_event) = self.queue.get_nowait()
                    # update the time of the scenario
                    self.t = round(time, self._time_rounding)
                    #self._counts_over_time[int(self.t)] = self._count_per_strains
                    # pass the event to the event handler
                    event_handler(n_event, get_neighbours)
                    # the new time is after the checking time
                    if self.t >= t_next_bin:
                        # check for the condition
                        for strain_id in surviving_strain_ids:
                            if not self.current_view.count(strain_id):
                                break
                except Empty:
                    self.log[round(self.t, self._log_time_rounding)] = copy(self.current_view)
                    break
            """
        # if we have a halt condition this part will conduct the simulation
        if with_halt_condition:
            halt = False
            # should the run be with explicit logs
            with_logging = params.get('explicit', False)  # defaults to False
            # work the event queue until we either hit t_stop or the halt condition
            while self.t < t_stop and not halt:
                try:
                    # get the next event
                    (time, n_event) = self.queue.get_nowait()
                    # update the time of the scenario
                    self.t = round(time, self._time_rounding)
                    #self._counts_over_time[int(self.t)] = self._count_per_strains
                    # pass the event to the event handler
                    event_handler(n_event, get_neighbours)
                    # the new time is after the checking time
                    if self.t >= t_next_bin:
                        if with_halt_condition:
<<<<<<< endemic/nw_spread/Spreading.py
                            self.log[round(self.t, self._log_time_rounding)] = copy(self.current_view)
                        while self.t >= t_next_bin:
                            t_next_bin += dt
=======
                            self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                        t_next_bin += dt
>>>>>>> endemic/nw_spread/Spreading.py
                        # check if we are in quasistable state (QSS) if yes, stop the sim
                        if self.quasistable(focus_strain_ids, surviving_strain_ids):
                            halt = True
                            # if we were not logging, write to the log now.
                            if not with_logging:
                                self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                # if no more events are to handle the sim is over (obviously)
                except Empty:
                    self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                    break
        # if we are in the case where a strain should build up its prevalence
        elif 'building_up' in params:
            # handle the stop condition
            strains_stop_cond = params.pop('building_up')
            if 'absolute' in strains_stop_cond:
                abs_cond = strains_stop_cond.pop('absolute')
            else:
                abs_cond = None
            if 'relative' in strains_stop_cond:
                rel_cond = strains_stop_cond.pop('relative')
            else:
                rel_cond = None
            if abs_cond:
                abs_id_stop_cond = {self.pathogen.ids[name]: abs_cond[name] for name in abs_cond}
                #building_status = {_id: self._count_per_strains[_id] for _id in abs_id_stop_cond}
                # check if the stop condition is already met
                if any(self._count_per_strains[_id] >= abs_id_stop_cond[_id] for _id in abs_id_stop_cond):
                    return 0
            if rel_cond:
                rel_id_stop_cond = {}
                for name in rel_cond:
                    _id = self.pathogen.ids[name]
                    if 'total' in rel_cond[name][1]:
                        if rel_cond[name][1] == 'total_infected':
                            # put the condition on the total set of present pathogens
                            rel_id_stop_cond[_id] = (rel_cond[name][0], self.pathogen.ids.values())
                            # if building up criterion is met already, stop the phase
                            if self._count_per_strains[_id] >= rel_id_stop_cond[_id][0] * sum(
                                    [self._count_per_strains[path_id] for path_id in self.pathogen.ids.values()]
                            ):
                                return 0
                        else:
                            rel_id_stop_cond[_id] = (rel_cond[name][0], 'all')
                            # if building up criterion is met already, stop the phase
                            if self._count_per_strains[_id] >= rel_id_stop_cond[_id][0] * self.contact_structure.n:
                                return 0
                    else:
                        rel_id_stop_cond[_id] = (
                            rel_cond[name][0], self.pathogen.ids[rel_cond[name][1]]
                        )
                        # if the condition is matched already end phase
                        if self._count_per_strains[
                            _id
                        ] >= rel_id_stop_cond[_id][0] * self._count_per_strains[rel_id_stop_cond[_id][1]]:
                            return 0
                            #rel_id_stop_cond = {
                            #    self.pathogen.ids[name]: (rel_cond[name][0], self.pathogen.ids[rel_cond[name][1]])
                            #    for name in rel_cond
                            #}
                            #if any(
                            #        self._count_per_strains[_id] >=
                            #        rel_id_stop_cond[_id][0] * self._count_per_strains[rel_id_stop_cond[_id][1]]
                            #        for _id in rel_id_stop_cond):
                            #    return 0
            # clarify which type of condition is active and define appropriate tests:
            # to do: revert back to using self._count_per_strains as soon as count_per_strains is reliable again
            if not abs_cond and rel_cond:
                def test_cond(self):
                    for s_id in rel_id_stop_cond:
                        condition = rel_id_stop_cond[s_id][1]
                        ref_val = condition if type(condition) in [float, int] else (
                            sum(
                                self.current_view.count(strain_id) for strain_id in self.pathogen.ids.values()
                                #self._count_per_strains[strain_id] for strain_id in self.pathogen.ids.values()
                            ) if type(condition) is list else self.contact_network.n
                        )
                        if self.current_view.count(s_id) >= rel_id_stop_cond[s_id][0] * ref_val:
                            #if self._count_per_strains[s_id] >= rel_id_stop_cond[s_id][0] * ref_val:
                            print 'reached fraction', self._count_per_strains, [
                                self.current_view.count(i) for i in self.pathogen.ids.values()
                            ]
                            # if we stop, we need to provide a new starting time for the next phase
                            self._phase_passon = {'t_start': self.t}
                            return 1
                    #if any(
                    #        self._count_per_strains[_id] >=
                    #        rel_id_stop_cond[_id][0] * self._count_per_strains[rel_id_stop_cond[_id][1]]
                    #        for _id in rel_id_stop_cond):
                    #    self._phase_passon = {'t_start': self.t}
                    #    return 1
                    else:
                        return 0
            elif not rel_cond and abs_cond:
                def test_cond(self):
                    if any(self._count_per_strains[_id] >= abs_id_stop_cond[_id] for _id in abs_id_stop_cond):
                        self._phase_passon = {'t_start': self.t}
                        return 1
                    else:
                        return 0
            else:
                def test_cond(self):
                    if any(
                                    self._count_per_strains[_id] >=
                                            rel_id_stop_cond[_id][0] * self._count_per_strains[rel_id_stop_cond[_id][1]]
                                    for _id in rel_id_stop_cond):
                        self._phase_passon = {'t_start': self.t}
                        return 1
                    if any(self._count_per_strains[_id] >= abs_id_stop_cond[_id] for _id in abs_id_stop_cond):
                        self._phase_passon = {'t_start': self.t}
                        return 1
                    else:
                        return 0
            # run the phase:
            if 'explicit' in params and params['explicit'] is True:
                while self.t < t_stop:
                    try:
                        (time, n_event) = self.queue.get_nowait()
                        self.t = round(time, self._time_rounding)
                        event_handler(n_event, get_neighbours)
                        if self.t >= t_next_bin:
<<<<<<< endemic/nw_spread/Spreading.py
                            self.log[round(self.t, self._log_time_rounding)] = copy(self.current_view)
                            while self.t >= t_next_bin:
                                t_next_bin += dt
=======
                            self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                            t_next_bin += dt
>>>>>>> endemic/nw_spread/Spreading.py
                        if test_cond(self):
                            return 0
                    except Empty:
                        self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                        break
            else:
                while self.t < t_stop:
                    try:
                        (time, n_event) = self.queue.get_nowait()

                        self.t = round(time, self._time_rounding)
                        event_handler(n_event, get_neighbours)
                        if test_cond(self):
                            return 0
                    except Empty:
                        self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                        break
        # if there was neither a halt condition nor a building_up, this part will conduct the simulation
        else:
            if 'explicit' in params and params['explicit'] is True:
                while self.t < t_stop:
                    try:
                        (time, n_event) = self.queue.get_nowait()
                        self.t = round(time, self._time_rounding)
                        event_handler(n_event, get_neighbours)
                        if self.t >= t_next_bin:
<<<<<<< endemic/nw_spread/Spreading.py
                            self.log[round(self.t, self._log_time_rounding)] = copy(self.current_view)
                            while self.t >= t_next_bin:
                                t_next_bin += dt
=======
                            self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                            t_next_bin += dt
>>>>>>> endemic/nw_spread/Spreading.py
                    except Empty:
                        self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                        break
            else:
                while self.t < t_stop:
                    try:
                        (time, n_event) = self.queue.get_nowait()
                        self.t = round(time, self._time_rounding)
                        event_handler(n_event, get_neighbours)
                    except Empty:
                        self.log[round(self.t, self._log_time_rounding)] = self.get_current_view
                        break
        #print 'treatment', with_treatment
        return 0

    def _parameter_alternation(self, alternations):
        """
        E.g. {wild_type: {'transmission_rate': 1}} will set the transmission rate of the wild_type to 1

        :param alternations: a dictionary holding stain names as keys and a dictionary as values.
            each of the value-dictionaries can contain rates as keys and the new desired value as values.
        :return:
        """
        for strain_name in alternations:
            # get its id
            its_id = self.pathogen.ids[strain_name]
            if 'transmission_rate' in alternations[strain_name]:
                new_rate = alternations[strain_name]['transmission_rate']
                self.pathogen.trans_rates[its_id] = new_rate
                self.pathogen.trans_dists[its_id] = Distro(new_rate ** (-1), 10000)
            if 'recover_rate' in alternations[strain_name]:
                new_rate = alternations[strain_name]['recover_rate']
                self.pathogen.rec_rates[its_id] = new_rate
                self.pathogen.rec_dists[its_id] = Distro(new_rate ** (-1), 10000)
                #to do: finish with recover_type
            if 'selection_rate' in alternations[strain_name]:
                new_rates = alternations[strain_name]['selection_rate']
                self.pathogen.update_selection(concerns=its_id, new_rates=new_rates)

    # to do: method should create list of instantaneous infection events then let the list be digested by appropriate
    # event handlers.
    def _sync_event_queue(self, mode='keep'):
        """
        This method will consume and recreate the event queue based on the status on self.current_view and the mode
        specified in the mode argument.


        :type mode: str; list; dict
        :param mode: either a string ('keep', 'reset') or a dict {'wild_type': 'keep'/'reset', ...} or a
            list [0, 1, ...] indicating the whether to reset (1) or to keep (0) the recover time for each node
        :return:
        """
        # get the reset mode
        general, pathogen_specific, node_specific = False, False, False
        if type(mode) is list:
            node_specific = True
        elif type(mode) is dict:
            pathogen_specific = True
        else:  # node_specific == True
            general = True
        event_queue = []
        while True:
            try:
                event = self.queue.get_nowait()
            except Empty:
                break
            else:
                # filter infection events
                if event[1][1] == -1:  # just take the recover events
                    event_queue.append(event)
        nodes_to_deal = map(lambda x: 1 if x != -1 else 0, self.current_view)
        while len(event_queue):
            event = event_queue.pop()
            # get the considered node
            node = event[1][0]
            if nodes_to_deal[node]:  # if the nodes new state is not susceptible:
                nodes_to_deal[node] = 0  # this node does not need further attention
                # get the nodes token
                token_id = self.current_view[node]
                # create recover time and infection events for the neighbours if the token is not -1 (susceptible)
                # get the nodes recover time
                if general:
                    if mode != 'keep':
                        # could we not just take it's old recover time?
                        # no, as the recover time will change for individuals that were treated before but not anymore
                        recover_time = self.pathogen.rec_dists[token_id].get_val()
                    else:
                        recover_time = event[0] - self.t
                elif pathogen_specific:
                    pathogen_name = self.pathogen.names[token_id]
                    if pathogen_name in mode and mode[pathogen_name] == 'reset':
                        recover_time = self.pathogen.rec_dists[token_id].get_val()
                    else:
                        recover_time = event[0] - self.t
                else:
                    if mode[node]:  # 1 means reset
                        recover_time = self.pathogen.rec_dists[token_id].get_val()
                    else:
                        recover_time = event[0] - self.t
                self.queue.put_nowait(Event(self.t + recover_time, node, -1, True))

                nn = n_copy(self.contact_structure.nn[node])  # get its nearest neighbours
                if nn.size:  # if he has some neighbours, get the times at which the get infected
                    inf_times = self.pathogen.trans_dists[token_id].v_get(nn)
                else:  # if he has no neighbours, make empty list
                    inf_times = array([])
                nn = nn[inf_times < recover_time]  # keep only the neighbours in nn which have an infection time
                # smaller than the nodes recover_time
                inf_times = inf_times[inf_times < recover_time]  # same thing for the infection times
                for x in xrange(inf_times.size):  # put all the infection events of neighbours into the queue
                    self.queue.put_nowait(Event(self.t + inf_times[x], nn[x], token_id, True, node))
        if nodes_to_deal.count(1):
            # ToDo: this is not a good strategy as it will break an SI model. Don't drop mutation events in the queue
            #   and change each mutation event to a recovery.
            while 1 in nodes_to_deal:
                node = nodes_to_deal.index(1)
                nodes_to_deal[node] = 0  # this node will be dealt with next time
                token_id = self.current_view[node]
                if token_id != -1:
                    recover_time = self.pathogen.rec_dists[token_id].get_val()
                    self.queue.put_nowait(Event(self.t + recover_time, node, -1, True))
                    nn = n_copy(self.contact_structure.nn[node])  # get its nearest neighbours
                    if nn.size:  # if he has some neighbours, get the times at which the get infected
                        inf_times = self.pathogen.trans_dists[token_id].v_get(nn)
                    else:  # if he has no neighbours, make empty list
                        inf_times = array([])
                    nn = nn[inf_times < recover_time]  # keep only the neighbours in nn which have an infection time
                    # smaller than the nodes recover_time
                    inf_times = inf_times[inf_times < recover_time]  # same thing for the infection times
                    for x in xrange(inf_times.size):  # put all the infection events of neighbours into the queue
                        self.queue.put_nowait(Event(self.t + inf_times[x], nn[x], token_id, True, node))

    # to do: this method needs some make over ( is not and should not be used at the moment )
    # - self._counts_over_time is not properly defined anymore
    def quasistable(self, quasi_stable_strain_ids=None, surviving_strain_ids=None):
        """
        Stability check.
        If stable return True, else return False
        """
        if quasi_stable_strain_ids is not None:
            i_1 = int(self.t / 3.)
            i_2 = 2 * i_1
            max_diff = n_max(absolute(
                divide(
                    n_sum(self._counts_over_time[i_1:i_2], axis=0), n_sum(self._counts_over_time[i_2:], axis=0)
                )[quasi_stable_strain_ids]
            ))
            if abs(1 - max_diff) >= 0.02:
                return False
            else:
                print 'quasistable at t= ', self.t
                return True
        if surviving_strain_ids is not None:
            if not count_nonzero(self._counts_over_time[int(self.t)][surviving_strain_ids]):
                print 'protected strain died out at t= ', self.t
                return True
            else:
                return False
        return False

    @property
    def get_outcome(self):
        """
        This function should be called at the end of each phase.
        It computes all the necessary properties and returns them.


        Should be strain specific:
        {time1:
            {'network': {'n': size, 'degree_count': {0: how many nodes with deg 0 }}
             'wild_type':
                {'count': total number of infected
                 'degree_count': {0:...,    will count the nbr of infected indivs per degree
                 'acquired': how many acquired the type through mutation
                 'degree_acquired': {0: how many node of deg 0 acquired through mutation
                }
             'mutant': ...
            },
         time2: ...
         }
        """

        _output = {
            'network': {
                'n': self.contact_structure.n,
                # 'degree_count': degree_count,
            }
        }
        if self.contact_structure.is_static:
            degrees = []
            for node in xrange(self.contact_structure.n):
                degrees.append(
                    len(
                        self.contact_structure.nn[node]
                    )
                )
            degree_count = {}
            observed_degrees = list(set(degrees))
            # Get nodes per degree
            nodes_per_degree = {deg: [] for deg in observed_degrees}
            for node_id in xrange(self.contact_structure.n):
                nodes_per_degree[degrees[node_id]].append(node_id)
            for a_degree in observed_degrees:
                degree_count[a_degree] = len(nodes_per_degree[a_degree])
            _output['network']['degree_count'] = degree_count
        else:
            # to do: what output for a dynamic network?
            pass
        # Run for each strain (we could also do the for strain_id ... inside the 2nd for loop below
        for strain_id in self.pathogen.names.keys():
            name = self.pathogen.names[strain_id]
            _output[name] = {}
            count = self.current_view.count(strain_id)
            _output[name]['count'] = copy(count)
            strain_acquired = 0
            if self.contact_structure.is_static:
                strain_degree_count = {}
                strain_degree_acquired = {}
                for a_degree in observed_degrees:
                    strain_degree_count[a_degree] = 0
                    strain_degree_acquired[a_degree] = 0
                    for node_id in nodes_per_degree[a_degree]:
                        if self.current_view[node_id] == strain_id:
                            strain_degree_count[a_degree] += 1
                            if self.current_infection_type[node_id] == 1:
                                strain_acquired += 1
                                strain_degree_acquired[a_degree] += 1
                _output[name]['degree_count'] = copy(strain_degree_count)
                _output[name]['acquired'] = copy(strain_acquired)
                _output[name]['degree_acquired'] = copy(strain_degree_acquired)
            else:
                for node_id in xrange(self.contact_structure.n):
                    if self.current_view[node_id] == strain_id:
                        if self.current_infection_type[node_id] == 1:
                            strain_acquired += 1
                _output[name]['acquired'] = copy(strain_acquired)
        return _output

    def _old_get_outcome(self):
        survived = {}
        degree_infected_fraction = {}
        acquire_type = {}
        degree_count = {}  # just a count for each degree
        degree_spec_acquire_count = {}  # Keys are the strain names, values are dict{degree: mutant count}
        degrees = []
        # degrees = {}
        for node in xrange(self.contact_structure.n):
            degrees.append(
                len(
                    self.contact_structure.nn[node]
                )
            )
            # degrees[node] = len(
            #     self.contact_structure.nn[node])  #to do: would make sense to define degree in contact_structure
        # observed_degrees = list(set(degrees.values()))
        observed_degrees = list(set(degrees))
        _degree_nodes = {deg: [] for deg in observed_degrees}  # for each degree a list of node ids
        infection_count_per_degree = {deg: {} for deg in observed_degrees}
        mutant_count_per_degree = {deg: {} for deg in observed_degrees}
        # the degree for each host
        # note: could also look at the average nearest neighbour degree...
        times = self.log.keys()  # to do: we need just the last time here, so just use self.t?!
        times.sort()
        # fraction of infected hosts
        survived[-1] = self.current_view.count(-1) / float(self.contact_structure.n)
        for strain_id in self.pathogen.ids.values():
            the_count = self.current_view.count(strain_id)
            survived[self.pathogen.names[strain_id]] = the_count / float(self.contact_structure.n)
        for node_id in xrange(self.contact_structure.n):  #to do: run through the strains rather than the nodes
            _degree_nodes[degrees[node_id]].append(node_id)
            state = copy(self.current_view[node_id])  #issue: probably not what causes the deviation...
            #state = last[node_id]
            the_strain = state
            if the_strain != -1:
                the_strain = self.pathogen.names[the_strain]
            try:
                infection_count_per_degree[degrees[node_id]][the_strain] += 1
            except KeyError:
                for strain_name in self.pathogen.names.values():
                    infection_count_per_degree[degrees[node_id]][strain_name] = 0.0
                infection_count_per_degree[degrees[node_id]][the_strain] = 1
        # Get degree specific acquire type
        degree_spec_acquire_type = []
        for degree in degrees:
            for a_node in _degree_nodes[degree]:
                degree_spec_acquire_type[degree].appennd(self.current_infection_type[a_node])
        degree_count = Counter(degrees)
        for degree in infection_count_per_degree:
            for strain_id in infection_count_per_degree[degree]:
                try:
                    degree_infected_fraction[degree][
                        strain_id
                    ] = infection_count_per_degree[degree][strain_id] / float(degree_count[degree]) \
                        if degree_count[degree] != 0. else 0.
                except KeyError:
                    degree_infected_fraction[
                        degree
                    ] = {strain_id: infection_count_per_degree[degree][strain_id] / float(degree_count[degree])} \
                        if degree_count[degree] != 0. else {strain_id: 0.}
        return {'abundance': survived, 'abundance_degree_specific': degree_infected_fraction}

    def _update_phase_in_sim_log(self, **params):
        self.simulation_log['phases'][self.t] = params
        self.simulation_log['phases'][self.t]['network'] = self.contact_structure.info
        self.simulation_log['phases'][self.t]['pathogen'] = self.pathogen.info
        if self.treatment is not None:
            self.simulation_log['phases'][self.t]['treatment'] = self.treatment.info
        self.simulation_log['phases'][self.t]['acquire_type'] = copy(self.current_infection_type)
        return 0

    # ToDo: WRITE A METHOD THAT RENDERS SELF.LOG AND SELF.SIMULATION_LOG MORE READABLE

    @property
    def get_current_view(self):
        if not self.contact_structure.has_dynamic_nodes or \
                (self.contact_structure.has_dynamic_nodes and self._ignore_dyn_nodes_in_log):
            return copy(self.current_view)   # default behaviour
        else:
            # Set all nodes that exceeded their lifespan (given by nodes_end) to a special state -2
            new_current_view = map(
                lambda i: self.current_view[i] if self.contact_structure.nodes_end[i] > self.t else -2,
                xrange(len(self.current_view)))
            return copy(new_current_view)

    # to transform the priority queue holding the upcoming events into a pickleabel list
    def __getstate__(self):
        d = dict(self.__dict__)
        queue = d.pop('queue')
        event_queue_list = []
        while True:
            try:
                event_queue_list.append(queue.get_nowait())
            except Empty:
                break
        d['event_queue_list'] = event_queue_list
        return d

    # to load the pickled event list back into a priority queue
    def __setstate__(self, d):
        if 'event_queue_list' in d:
            event_queue_list = d.pop('event_queue_list')
            d['queue'] = PriorityQueue()
            while len(event_queue_list):
                d['queue'].put_nowait(event_queue_list.pop())
        self.__dict__.update(d)


# issue: Make this a static method of Spreading?
def Event(time, node, token, inf_event, source=None):
    """
    Arguments:
        - time: float, used to order the events in the priority queue
        - node: name of the affected host
        - token: The putative new infectious status of the host
        - inf_event: Whether the Event is an infection (True) or a mutation(False)
        - source: The id of the node where the infection came from
    """
    return time, (node, token, inf_event, source)