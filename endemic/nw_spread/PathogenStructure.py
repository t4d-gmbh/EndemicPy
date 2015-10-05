from RateDistribution import Distro

#this is passed as the scale parameter if rate == 0 - see Distro.fillup() for details.
NO_RATE_FLAG = None


class Pathogen():
    def __init__(self, strains):
        """
        A pathogen object can hold a collection of strains (object form Strain class).
        At some point it might even be possible that a pathogen object
            can create new _strains. But for now it just holds a collection
            of previously defined Strain objects
        
        Arguments:
            - strains: a list of Strain objects
        """
        self._strains = []  #list of ids
        self.n = 0  #gives the number of strains
        self.names = {}
        self.ids = {v: k for k, v in self.names.iteritems()}
        self.trans_rates = []
        self.trans_dists = []
        self.rec_rates = []
        self.rec_dists = []
        self.rec_types = []
        self.select_rates = []
        self._select_rates = []
        self.select_dists = {}  #NOTE: this is going to be rate specific
        for a_strain in strains:
            self._append_strain(a_strain)
        self._check_integrity()
        self._resolve_selection()
        self._get_selection_dists()

    def _append_strain(self, a_strain):
        its_id = len(self._strains)
        self.names[its_id] = a_strain.name
        self.ids[a_strain.name] = its_id
        self._strains.append(its_id)
        self.n = len(self._strains)
        self.trans_rates.append(a_strain.b)
        self.trans_dists.append(a_strain.b_dist)
        self.rec_rates.append(a_strain.g)
        self.rec_dists.append(a_strain.g_dist)
        self.rec_types.append(a_strain.recover_type)
        self._select_rates.append(a_strain.s)
        return 0

    #issue: this is to modify an existing strain = not used for now
    def _modify_strain(self, strain_id, new_strain):
        its_id = strain_id
        self.names[its_id] = new_strain.name
        self.ids[new_strain.name] = its_id
        self._strains.append(its_id)
        #self.n = len(self._strains)
        self.trans_rates[its_id] = new_strain.b
        self.trans_dists[its_id] = new_strain.b_dist
        self.rec_rates[its_id] = new_strain.g
        self.rec_dists[its_id] = new_strain.g_dist
        self.rec_types[its_id] = new_strain.recover_type
        self._select_rates[its_id] = new_strain.s
        #to do: selection still needs to be handled: _resolve_selection and _get_selection_dists
        return 0

    @property
    def info(self):
        return {
            'pathogens': self.ids.keys(),
            'rates': {
                'transmission': [self.trans_rates[self.ids[_name]] for _name in self.ids.keys()],
                'recover': [self.rec_rates[self.ids[_name]] for _name in self.ids.keys()],
                'selection': [self._select_rates[self.ids[_name]] for _name in self.ids.keys()]
            }
        }

    class MissingStrainError(Exception):
        pass

    def _check_integrity(self, ):
        """
        Check the integrity of the ensemble of _strains.
        """
        if len(list(set(self._strains))) != len(self._strains):
            raise self.MissingStrainError('Not all strains have distinct name entries. Please provide a set of strains\
             with all unique names')
        else:
            pass
        return 0

    def _resolve_selection(self):
        """
        This populates the self.select_rates list.

        self.select_rates gives for each strain (id is index) the list of mutation rates towards the other strains.
        :return:
        """
        names = self.ids.keys()
        for _id in self._strains:
            if 'Default' in self._select_rates[_id]:
                def_val = self._select_rates[_id].pop('Default')
            else:
                def_val = 0
            strain_array = [def_val for _ in xrange(self.n)]
            for name in self._select_rates[_id]:
                if name not in names:
                    raise self.MissingStrainError(
                        "In the selection rates of strain {0:s} appears a unknown strain: {1:s}".format(self.names[_id],
                                                                                                        name)
                    )
                strain_array[
                    self.ids[name]
                ] = round(self._select_rates[_id][name], 6)  #Note: The selection rates are rounded on 6 digits
            self.select_rates.append(strain_array)
        return 0

    def update_selection(self, concerns, new_rates):
        """
        This updates the self.select_rates list.

        self.select_rates gives for each strain (id is index) the list of mutation rates towards the other strains.
        :return:
        """
        self._select_rates[concerns] = new_rates
        names = self.ids.keys()
        for _id in self._strains:
            if 'Default' in self._select_rates[_id]:
                def_val = self._select_rates[_id].pop('Default')
            else:
                def_val = 0
            strain_array = [def_val for _ in xrange(self.n)]
            for name in self._select_rates[_id]:
                if name not in names:
                    raise self.MissingStrainError(
                        "In the selection rates of strain {0:s} appears a unknown strain: {1:s}".format(self.names[_id],
                                                                                                        name)
                    )
                strain_array[
                    self.ids[name]
                ] = round(self._select_rates[_id][name], 4)  #Note: The selection rates are rounded on 4 digits
            self.select_rates[_id] = strain_array
            self._get_selection_dists()
        return 0

    def _get_selection_dists(self):
        """
        This creates the self.selection_dists list.

        self.selection_dists holds for each selection rate existing in the scenario a Distro object to draw a seleciton
            time from.
        :return:
        """
        selection_rates = []
        for indiv_selections in self.select_rates:
            selection_rates.extend(indiv_selections)
        selection_rates = list(set(selection_rates))
        selection_rates.sort()
        for a_rate in selection_rates:
            if a_rate not in self.select_dists:
                try:
                    scale = a_rate ** (-1)
                except ZeroDivisionError:
                    scale = NO_RATE_FLAG  #set the flag for a rate of 0.
                self.select_dists[a_rate] = Distro(scale=scale, pre=10000)
        return 0

    @staticmethod
    def index_min(values):
        return min(xrange(len(values)), key=values.__getitem__)

    #to do: delete this. its the old version of the method below
    def _get_selected(self, strain_id, boosts=None):
        """

        :param strain_id: The id of the strain that is introduced
        :param boosts:
        :return:
        """
        potential_selection = self.select_rates[strain_id]  #list: index(strain_id), value(selection_rate)
        if boosts is None:
            selection_times = map(lambda x: self.select_dists[x].get_val(), potential_selection)
        else:
            #issue: look closer at how the boosts influence the selection time
            selection_times = map(lambda x: self.select_dists[x].get_val() * boosts[strain_id], potential_selection)
        i = self.index_min(selection_times)
        return i, selection_times[i]  #strain_id and time of selection

    #issue: the boost as well as the delay is not needed.
    def get_selected(self, strain_id):
        """

        :param strain_id: The id of the strain that is introduced
        :return:
        """
        potential_selection = self.select_rates[strain_id]  #list: index(strain_id), value(selection_rate)
        selection_times = map(lambda x: self.select_dists[x].get_val(), potential_selection)
        i = self.index_min(selection_times)
        return i, selection_times  #strain_id and time of selection


class Strain():
    def __init__(self, name, transmission_rate, recover_rate, recover_type=1, selection_rate=None, ):
        """
        This class defines a single strain.
        
        Arguments:
            - name: An unique identifier for a strain.
                NOTE: It is assumed that name=0 is the wild type.
            - transmission_rate: the rate of transmission for this strain
            - recover_rate: the rate of recovery from this strain.
            - recover_type: Gives the status of a host after recovery form this
                strain.
                Possible choices are:
                    - 0: The host becomes resistant to this strain -> SIR-model
                    - 1: The host becomes susceptible again -> SIS-model
                    - Default=1: SIS-model
                    Note: Here we can also implement partial resistance (0<recover_type<1)
            - selection_rate: The rate at which the strain mutates and gets selected.
                Default=None: The strain does not mutate
                The rate can either be given as a float or a dict.
                If it is a float, the strain mutates to each one of the active 
                    strains with equal probability and a rate given by the value.
                It it is a dict, the keys must correspond to a strain name and the
                    respective value gives the mutation rate towards this strain.
                    If 'Default' is given in the dict, all missing strains get the value of 'Default'.
                    E.g. selection_rate={'wild_type':0.5,...}
                    A selection rate of 0 means no selection.
        """
        self.name = name
        self.b = transmission_rate
        try:
            self.b_dist = Distro(self.b ** (-1), 10000)
        except ZeroDivisionError:
            self.b_dist = Distro(NO_RATE_FLAG, 10000)
        self.g = recover_rate
        try:
            self.g_dist = Distro(self.g ** (-1), 10000)
        except ZeroDivisionError:
            self.g_dist = Distro(NO_RATE_FLAG, 10000)
        if selection_rate:
            if type(selection_rate) is not dict:
                self.s = {'Default': selection_rate}
            else:
                self.s = selection_rate
        else:
            self.s = {'Default': 0}
        self.recover_type = recover_type

    #write get_val functions for the different rates...is presumably faster than calling the fct form the Distro class