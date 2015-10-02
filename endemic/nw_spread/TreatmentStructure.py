class Treatment():
    def __init__(self, therapies):
        """
        Note: Each pathogen strain can only be treated by 1 Therapy at once.
        
        Arguments:
            - _therapies is a list of Therapy object
        """
        self.therapies = therapies
        self.n = len(self.therapies)
        self.names = {}
        self.ids = {}
        #could call _append_therapy here rest is not needed
        for a_therapy in therapies:
            self._append_therapy(a_therapy)
        self._check_integrity()

    def _append_therapy(self, a_therapy):
        its_id = len(self.names)
        self.names[its_id] = a_therapy.name
        self.ids[a_therapy.name] = its_id
        return 0

    def add_new_treatment(self, treatment):
        """
        Method to add therapies from another treatment object
        :param treatment:
        :return:
        """
        for therapy in treatment.therapies:
            self.therapies.append(therapy)
            self._append_therapy(therapy)
        self.n = len(self.therapies)
        return 0

    @property
    def info(self):
        return {
            'names': [a_therapy.name for a_therapy in self.therapies],
            'factors': {
                'transmission': [a_therapy.drug.transmission_factor for a_therapy in self.therapies],
                'recovery': [a_therapy.drug.recover_factor for a_therapy in self.therapies],
                'selection': [a_therapy.drug.selection_factor for a_therapy in self.therapies]
            }
        }

    class MissingTherapyError(Exception):
        pass

    def _check_integrity(self, ):
        """
        Check the integrity of the ensemble of _strains.
        """
        if len(list(set(self.ids.keys()))) != len(self.ids.keys()):
            raise self.MissingTherapyError('Not all names have distinct name entries. Please provide a set of\
             names with all unique names')
        else:
            pass
        return 0


class Therapy():
    def __init__(self, name, drug, delay=0.0, treatment_proba=0.0):
        """
        This class defines different strategies for a treatment.
        
        Arguments:
            - name: A name for the therapy (each therapy name must be unique)
            - drug, an object from the drug class, the drug that is used for the
                treatment.
            - delay, Default=0.0 the time delay between infection and treatment
            - treatment_proba, Default=0.0, the probability of treating an infected host.
                it can either be a float (uniform probability) a list (proba for each node) or
                one of the following strings: ...
        """
        self.name = name
        self.delay = delay
        self.treatment_proba = treatment_proba
        self.drug = drug


class Drug():
    class InvalidSelectionBoostError(Exception):
        pass

    def __init__(self, name, transmission_factor=None, recover_factor=None, selection_factor=None):
        """
        Arguments:
            - name: A unique identifier for the drug
            - transmission_factor: Either a float or a dict giving the inhibitory factor(s) by
                which the transition rate of the targeted strain is multiplied.
                If a dict is provided, the key must be the name of a strain and the
                    value the alternation factor. If not all strain IDs are present
                    the absent strains are not affected by the drug.
            - recovery_factor: Ether a float or a dict giving the boost factor(s) by
                which the recovery rate of the treated host is multiplied.
            - selection_factor: Either a float or a dict giving the boost factors for mutation+selection.
                If a float is given, then all possible selection processes will get this boost.
                If a dict is given, (key: strain name, value the boost factor) then each strain will get its specific
                    boost factor. If 'Default' is given, then all missing strains get the default value. If 'Default'
                    is not given but strains are missing, then the default value is 1, meaning no change.
                    Note: the factor is considered independent from the source strain, i.e. a factor of 2 for 'mutant_2'
                        will double the mutation rate from any other strain to 'mutant_2'.

        """
        self.name = name
        #self.transmission_inhibitor = 1
        self.transmission_factor = {'Default': 1.}
        if transmission_factor is not None:
            if type(transmission_factor) is dict:
                for strain_name in transmission_factor:
                    self.transmission_factor[strain_name] = transmission_factor[strain_name]
            else:
                #self.transmission_inhibitor = transmission_inhibitor
                self.transmission_factor['Default'] = transmission_factor
        #self.recover_boost = 1
        self.recover_factor = {'Default': 1.}
        if recover_factor is not None:
            if type(recover_factor) is dict:
                #self.recover_boost = recover_boost
                for strain_name in recover_factor:
                    self.recover_factor[strain_name] = recover_factor[strain_name]
            else:
                #self.recover_boost = recover_boost
                self.recover_factor = {'Default': recover_factor}
        self.selection_boost = {'Default': 1.}
        if selection_factor is not None:
            if type(selection_factor) is int or type(selection_factor) is float:
                self.selection_boost['Default'] = selection_factor
            elif type(selection_factor) is dict:
                for strain_name in selection_factor.keys():
                    self.selection_boost[strain_name] = selection_factor[strain_name]
            else:
                raise self.InvalidSelectionBoostError('selection_boost is of the wrong format. See docstring for info.')
        self.selection_factor = self.selection_boost  #issue: this is ugly, self.selection_boost is not needed
        #self.selection_factor = {}
        #for name in self.selection_boost:
        #    self.selection_factor[name] = self.selection_boost[name] ** (-1)