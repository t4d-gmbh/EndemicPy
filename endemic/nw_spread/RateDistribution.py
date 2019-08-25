from numpy import vectorize, random, apply_along_axis, array, float64
from Queue import Empty
from Queue import Queue as SimpleQueue

#this is the time it takes if a rate of 0 is given
MAX_LIM = 100000


def no_mut(dummy, length):
    return array([MAX_LIM] * length, dtype=float64)


class Distro(object):
    def __init__(self, scale, pre=10):
        """
        This class holds a queue of times drawn from an exponential 
            distribution with a specified scale.
        
        Arguments:
        
            - scale: The scale parameter for the exponential distribution.
            - pre: Predefined size of the queue. Default=10
        """
        self.scale = scale
        self.pre = pre
        self.queue = SimpleQueue(maxsize=pre + 1)
        self.v_put = vectorize(self.queue.put_nowait)
        #the exponential dist is not defined for a rate of 0
        #therefore if the rate is 0 (scale is None then) huge times are set
        if self.scale in [None, 0]:
            self.scale = 0
            self.draw_fct = no_mut
        else:
            self.draw_fct = random.exponential
        #fillup the queue
        self.fillup()
        # there was: (new version compatible with pickeling see method below)
        self.v_get = vectorize(self.get_val)

    def fillup(self):
        self.v_put(self.draw_fct(self.scale, self.pre))
        return 0

    def get_val(self, a=None):
        """
        Function returning a value drawn form the exponential distribution.
        """
        try:
            return self.queue.get_nowait()
        except Empty:
            self.fillup()
            return self.queue.get_nowait()
        
    #def v_get(self, an_array):
    #    #return map(self.get_val, xrange(an_array.size))
    #    return apply_along_axis(self.get_val, 0, an_array)
    # to transform the priority queue holding the upcoming events into a pickle-abel list
    
    def __getstate__(self):
        d = dict(self.__dict__)
        queue = d.pop('queue')
        # v_put is not pickle-able
        del d['v_put']
        del d['v_get']
        simple_queue_list = []
        while True:
            try:
                simple_queue_list.append(queue.get_nowait())
            except Empty:
                break
        d['simple_queue_list'] = simple_queue_list
        return d

    # to load the pickled event list back into a priority queue
    # TO DO: oh boy this is some horrible patch
    def __setstate__(self, d):
        if 'simple_queue_list' in d:
            event_queue_list = d.pop('simple_queue_list')
            d['queue'] = SimpleQueue(maxsize=d['pre'] + 1)
            while len(event_queue_list):
                d['queue'].put_nowait(event_queue_list.pop())
        self.__dict__.update(d)
        self.__dict__['v_put'] = vectorize(self.queue.put_nowait)
        #d['v_put'] = vectorize(d['queue'].put_nowait)
        #self.__dict__.update(d)
        self.__dict__['v_get'] = vectorize(self.get_val)
        if self.scale is None:
            self.scale = 0
            self.queue = SimpleQueue(maxsize=self.pre + 1)
            self.v_put = vectorize(self.queue.put_nowait)  # this is specific to the queue, thus reinit here
            self.draw_fct = no_mut
            self.fillup()
