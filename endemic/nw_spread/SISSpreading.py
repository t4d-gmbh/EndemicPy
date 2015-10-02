#To Do: This is old and should be removed.

# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import numpy.random as n_random
from Queue import Empty
from Queue import Queue as simple_queue
from Queue import PriorityQueue as p_queue

# <codecell>

class SIS_spreading(object):
    def __init__(self, graph, g, b, spreading_order = 'advance', **params):
        """graph is an igraph.Graph object
           g, the recovery rate, is a dictionary, giving the rate for each viral strain
           b, the transmission rate, dito a dict. giving the rate for each viral strain
           spreading_order is a string, indicating the type of spreading. Defautl is 'advance'.
               possible choices: 'advance': a new strain emerges once the previous reached endemic susceptible
                                 'parallel': all _strains start together
                                 'sequential': ...
           params are additional parameters for the simulation:
               t_max: gives the maximal time after which the simulations stops. Default is 1000
               t_burn: gives the time window during which a virus can initially spread.
         """
        self.G = graph
        self.spreading_order = spreading_order
        #initialize population
        self.status = np.empty(self.G.N, dtype = np.int)
        self.status.fill(0)
        #initialize other parameters
        self.parameters = params
        self._update_parameters(g, b)
        #init the 'time'
        self.t = 0
        #inti the result dict
        self.prevalence = {}
        #form {time: [count_0,count_1,...]} with count_0 the healty, count_1 infected by strain 1...
        
    def _update_parameters(self, g, b):
        #gives the number of viral _strains
        self.strains = max(len(g.keys()), len(b.keys()))
        self.parameters['g'] = g
        self.parameters['b']= b
        if 'dt' not in self.parameters:
            try:
                self.dt = 0.1*g[1]
            except KeyError:
                self.dt = 1
        else:
            self.dt = self.parameters['dt']
        if 't_burn' not in self.parameters:
            try:
                self.t_burn = 100/float(g[1])
            except KeyError:
                self.t_burn = 100
            self.parameters['t_burn'] = self.t_burn
        if 't_max' not in self.parameters:
            self.t_max = 1000
            self.parameters['t_max'] = self.t_max
        #variable that holds the list of active _strains
        #self.active_strains = np.array(range(1,self._strains+1))
        self.active_strains = np.arange(1,self.strains+1)

    def _starting_event(self, strain = 1, potential_sources = False):
        """choses an initial note to start the spreading
           strain: gives the id for the strain
           potential_sources: list of all potential souces (node ids)
        """
        if not potential_sources:
            start_index = np.random.randint(0,self.G.N)
        else:
            start_index = potential_sources[np.random.randint(0,len(potential_sources))]
        self.status[start_index] = strain
        return 0
        
    def reset_infection(self,strains_to_reset = False):
        """resets one or several _strains
        """
        if not strains_to_reset:
            strains_to_reset = np.arange(1,self.strains+1)
        self.status = np.where(np.in1d(self.status,strains_to_reset),0,self.status)
        self.t = 0
        self.t_burn = self.parameters['t_burn']
        self.t_max = self.parameters['t_max']
        self.prevalence = {}
        return 0
    

# <codecell>

class SIS_continuous(SIS_spreading):
    """Event based implementation of the SIS process
        graph is an igraph.Graph
        g,b are the recovery, respectively transmission drates in form of dictionaries:
            Example: g = {1:rate_1,2:rate_2, ...}
        the number of items in g and b should match, they determine the number of different virus _strains.
        spreading_order determines the structure of the spreading process:
            advance (default): 
        params are optional additional parameters:
            t_burn: burn in period (default 100/g)
            t_max: maximal simulation time (default t_max = 1000)
        """
    def __init__(self, graph, g, b, spreading_order = 'advance', **params):
        """g and b are of the form: g = {1: scale, 2: scale}"""
        SIS_spreading.__init__(self,graph,g,b,spreading_order,**params)
        #the current time in the simulation
        self.t = 0
        #create the priorityqueue
        self.queue = p_queue()
        #vectorized versions of the get and put command
        self.queue.v_get = np.vectorize(self.queue.get_nowait)
        self.queue.v_put = np.vectorize(self.queue.put_nowait)
        self._update_rates(g,b)
            
    def _update_rates(self,g,b):
        self.g = {}
        self.b = {}
        for key in self.parameters['g']:
            try:
                g_scale = g[key]**(-1)
            except ZeroDivisionError:
                g_scale = 0
            self.g[key] = distro(g_scale)
        for key in self.parameters['b']:
            try:
                b_scale = b[key]**(-1)
            except ZeroDivisionError:
                b_scale = 0
            self.b[key] = distro(b_scale)
    
    def starting_event(self, strain = 1, potential_sources = False):
        """choses an initial note to start the spreading
           strain: gives the id for the strain
           potential_sources: list of all potential souces (node ids)
        """
        if not potential_sources:
            self._starting_event(strain = strain)
        else:
            self._starting_event(strain = strain, potential_sources = potential_sources)
        self._queue_init(strain)
        return 0
    
    def _queue_init(self, focus_strain = 1):
        #find the initially changed node.
        node = np.where(self.status == focus_strain)[0][0]
        strain = focus_strain
        recover_time = self.g[strain].get_val()
        #check if the node is infected with something else already
        #if so, reset the recover Event.
        if self.queue.qsize():
            found = False
            put_back_events = []
            while not found:
                an_event = self.queue.get_nowait()
                if an_event[1][0] != node:
                    put_back_events.append(an_event)
                else:
                    found = True
            map(lambda x:self.queue.put_nowait(x),put_back_events)
        self.queue.put_nowait(event(self.t + recover_time,node,0))
        #get neighbors indexes
        #nn = np.where(self.g_dist.nn[node] == 1)[0]
        nn = np.array(self.G.NN[node])
        #if strain == 2:
            #print 'got ', nn.size
        if nn.size:
            #get an array of recovery times
            inf_times = self.b[strain].v_get(nn)
            #filter for the events that occurre at time inf_time < recover_time
            nn = nn[inf_times<=recover_time]
            #print 'will inf ', nn.size
            #print nn
            #print [self.status[nnn] for nnn in nn]
            inf_times = inf_times[inf_times<=recover_time]
            map(lambda x:self.queue.put_nowait(event(self.t +inf_times[x],nn[x],strain)),xrange(int(inf_times.size)))
            return 0
        return 0
            
    def handle_event(self,an_event):
        (node,strain) = an_event
        if not strain: #and not diseas_immunity[node]:
            self.status[node] = 0
        #in case it is an infection Event
        else:
            #if the node can be infected
            if not self.status[node]:
                #create the recover Event
                recover_time = self.g[strain].get_val()
                self.queue.put_nowait(event(self.t + recover_time, node, 0))
                self.status[node] = strain
                #nn = np.where(self.g_dist.nn[node] == 1)[0]
                #if len(self.g_dist.nn[node]):
                nn = np.array(self.G.NN[node])
                #get an array of recovery times
                inf_times = self.b[strain].v_get(nn)
                #filter for the events that occurre at time inf_time < recover_time
                #print 'has ', nn.size, 'neighbours'
                nn = nn[inf_times<=recover_time]
                inf_times = inf_times[inf_times<=recover_time]
                #diseas_immunity[nn] +=inf_times
                #for each nn put a recovery time, a disease
                for x in xrange(int(inf_times.size)):
                    self.queue.put_nowait(event(self.t + inf_times[x],nn[x],strain))
                #map(lambda x:self.queue.put_nowait(Event(inf_times[x],nn[x],strain)),range(int(inf_times.size)))
            else:
                pass
        return 0
            
    def run(self, stop_time, t_burn, focus_strain = [1], stop_cond = True, log = True):
        state = False
        last_time = 0
        counter = 0
        start_time = self.t
        #print 'start time', start_time
        #diseas_immunity = np.empty()
        #diseas_immunity.fill(0)
        while self.t < stop_time:
            try:
                (time,n_event) = self.queue.get_nowait()
                #diseas_immunity - (time-self.t)
                #diseas_immunity = np.where(diseas_immunity >=0,diseas_immunity,0)
                self.t = time
                self.handle_event(n_event)
                if log:
                    if self.t - last_time >= self.dt:
                        counter += 1
                        self.prevalence[self.t] = list(np.bincount(self.status, minlength=self.strains+1))
                        last_time = self.t
                        if counter >= 20:
                            #print self.status
                            #print 't',self.t, 't_burn', t_burn
                            counter = 0
                            if any([self.prevalence[self.t][fs] == 0 for fs in focus_strain]):
                                state = 'extinct'
                                raise Empty
                            #if self.status.sum == 0:
                            #    break
                            if stop_cond and self.t >= t_burn:
                                #print 'check stop cond'
                                one3th_t = start_time + (self.t-start_time)*0.333
                                two3th_t = start_time + (self.t-start_time)*0.666
                                #print 'first', one3th_t
                                #print 'second', two3th_t
                                #print 't', self.t
                                items = self.prevalence.items()
                                #print map(lambda x:x[1],filter(lambda x:x[0],items))
                                last_avgs = np.average(map(lambda x:x[1],filter(lambda x:x[0]>=two3th_t,items)),axis = 0)
                                middle_avgs = np.average(map(lambda x:x[1],filter(lambda x:x[0]>=one3th_t and x[0]<two3th_t, items)),axis =0)
                                for strain in focus_strain:
                                    #print last_avgs
                                    #print middle_avgs
                                    if abs(last_avgs[strain]/middle_avgs[strain] - 1)<= 0.02:
                                        state = 'endemic_' + str(strain)
                                        raise Empty
            except Empty:
                if not state:
                    state = 'extinct'
                print state
                break
        if not state:
            state = 'at_max_t'
        #print 'END', self.t, 'burn', self.t_burn
        self.prevalence[self.t] = list(np.bincount(self.status, minlength=self.strains+1))
        return state
      
    def new_trial(self, g, b, **distribution):
        self.reset_infection()
        self.queue = p_queue()
        self.G.new_realization(**distribution)
        self._update_parameters(g,b)
        self._update_rates(g,b)
        return 0
        
    def simple_scenario(self,g=None,b=None,stop_cond = 'to_endemic',focus_strain = [1],**distribution):
        """run a symulation with a single strain"""
        if not g:
            g = self.parameters['g']
        if not b:
            b = self.parameters['b']
        self.new_trial(g=g,b=b,**distribution)
        self.starting_event(strain = 1, potential_sources = False)
        if stop_cond != 'to_endemic':
            t_stop = delai
        else:
            t_stop = self.t_max
        final_state = self.run(t_stop,t_burn = self.t_burn, focus_strain = focus_strain)
        status = list(self.status)
        count_1 = [status.count(i) for i in xrange(max(status)+1)]
        deg_count_1 = [[0 for i in xrange(max(status)+1)] for d in xrange(max(self.G.degrees)+1)]
        for i in xrange(self.G.N):
            deg_count_1[self.G.degrees[i]][status[i]] += 1
        return (self.t,final_state,count_1,deg_count_1)

    def double_scenario(self,g=None,b=None,delai = 'to_endemic',**distribution):
        """run a symulation with several _strains"""
        if not g:
            g = self.parameters['g']
        if not b:
            b = self.parameters['b']
        ok = False
        while not ok:
            try:
                self.new_trial(g=g,b=b,**distribution)
                ok = True
            except ValueError:
                ok = False
        self.starting_event(strain = 1, potential_sources = False)
        if delai != 'to_endemic':
            t_stop = delai
        else:
            t_stop = self.t_max
        first_state = self.run(t_stop,t_burn = self.t_burn, focus_strain = [1])
        first_status = list(self.status)
        count_1 = [first_status.count(i) for i in xrange(max(first_status)+1)]
        deg_count_1 = [[0 for i in xrange(max(first_status)+1)] for d in xrange(max(self.G.degrees)+1)]
        for i in xrange(self.G.N):
            deg_count_1[self.G.degrees[i]][first_status[i]] += 1
        t_1 = self.t
        if 'endemic' in first_state or first_state == 'at_max_t':
            self.starting_event(strain = 2, potential_sources = filter(lambda x:self.status[x] == 1,range(self.G.N)))
            #print self.status
            second_state = self.run(self.t_max+self.t, t_burn = self.t_burn + self.t,focus_strain = [2])
            second_status = list(self.status)
            count_2 = [second_status.count(i) for i in xrange(max(second_status)+1)]
            deg_count_2 = [[0 for i in xrange(max(second_status)+1)] for d in xrange(max(self.G.degrees)+1)]
            for i in xrange(self.G.N):
                deg_count_2[self.G.degrees[i]][second_status[i]] += 1
            #return {'1st_state':first_state,'1st_status':first_status,'2nd_state': second_state, '2nd_status':second_status}
            return ((t_1,first_state,count_1,deg_count_1),(self.t,second_state,count_2,deg_count_2))
        else:
            print 'The first strain died, no second wave possible'
            return ((t_1,first_state,count_1,deg_count_1),(None,None,None,None))

# <markdowncell>

# ---
# 
# 
# 
# ---
# 
# 
# 
# ###For the Event based simulation:

# <codecell>

def event(time,node,status):
    return (time,(node,status))

# <codecell>

class distro(object): 
    def __init__(self, scale, pre = 10):
        self.scale = scale
        self.pre = pre
        self.queue = simple_queue(maxsize = pre+1)
        self.v_put = np.vectorize(self.queue.put_nowait)
        #the exponential dist is not defined for scale = 0
        #therefore set scale to an infinitesimal value if 
        #it is 0.
        if not self.scale:
            self.scale = 0.00001
        #fillup the queue
        self.fillup()
        self.v_get = np.vectorize(self.get_val)
        
    def fillup(self):
        self.v_put(n_random.exponential(self.scale,self.pre))
            
    def get_val(self,a = None):
        try:
            return self.queue.get_nowait()
        except Empty:
            self.fillup()
            return self.queue.get_nowait()

# <markdowncell>

# ---
# 
# General SIS model

# <codecell>

class multiSpread(object):
    """This will be the general class for SIS/SIR spreading of multiple _strains"""

# <markdowncell>

# ---
# 
# Strain model

# <codecell>

class Strain(object):
    """This class defines various possible _strains"""

# <markdowncell>

# ---
# Host model

# <codecell>

class Host(object):
    """Defines the _hosts on which a pathogen spreads."""
    #might be userfull to create a minimal class
    #Attributes:
    # - list of neighbours
    # - infection status
    # - treatment status
    # - resistant status
