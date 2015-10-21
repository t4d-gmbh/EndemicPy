from numpy import vectorize, array, float64
from scipy.stats import expon, norm, lomax
from Queue import Empty
from Queue import Queue as SimpleQueue

#this is the time it takes if a rate of 0 is given
MAX_LIM = 10000


def no_mut(loc=None, scale=None, size=10):
    return array([MAX_LIM] * size, dtype=float64)


class Distro(object):
    """
    This class holds a queue of random times drawn from a given distribution with a specified scale.
    
    Arguments:
    
        - scale: The scale parameter for the exponential distribution. (default=1)
        - pre: Predefined size of the queue. (default=10)
        - loc: location of the distribution (default=0)
        - alpha: power law exponent (default=1), only for the lomax distribution
        - distribution type: string of the name of a distribution from scipy.stats
            "exponential" (default) with parameter scale:
                pdf(x, scale) = (1/scale) * exp(-x/scale)
                
            "gaussian":
                pdf(x, scale, loc) = |(1/scale) * norm((x-loc)/scale)|
                note : returns the absolute value
                
            "lomax" power-law distribution:
                pdf(x, alpha) = alpha / (1+x)**(alpha+1)
                for x >= 0, alpha > 0
                the minimum value is controlled by `loc`
            
    """
    def __init__(self, scale=1, pre=10, loc=0, alpha=1, distribution_type="exponential"):
       
        self.scale = scale
        self.loc = loc
        self.pre = pre
        self.alpha = alpha
        self.queue = SimpleQueue(maxsize=pre + 1)
        self.v_put = vectorize(self.queue.put_nowait)
        self.distribution_type = distribution_type
        
        
        #the exponential dist is not defined for a rate of 0
        #therefore if the rate is 0 (scale is None then) huge times are set
        if self.scale is None or self.scale == 0:
            self.draw_fct = no_mut
        else:
            if distribution_type == "exponential":                
                self.draw_fct = expon.rvs
            elif distribution_type == "gaussian":
                self.draw_fct = norm.rvs
            elif distribution_type == "lomax":
                self.draw_fct = lomax.rvs
            else:
                raise self.DistributionTypeError('This distribution type is not implemented yet.')
                
        #fillup the queue
        self.fillup()
        # there was: (new version compatible with pickeling see method below)
        self.v_get = vectorize(self.get_val)

    class DistributionTypeError(Exception):
        pass
    
    def fillup(self):
        if self.distribution_type == "lomax":
            self.v_put(abs(self.draw_fct(self.alpha, loc=self.loc, scale=self.scale, size=self.pre)))
        else:
            self.v_put(abs(self.draw_fct(loc=self.loc, scale=self.scale, size=self.pre)))
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
    def __setstate__(self, d):
        if 'simple_queue_list' in d:
            event_queue_list = d.pop('simple_queue_list')
            d['queue'] = SimpleQueue()
            while len(event_queue_list):
                d['queue'].put_nowait(event_queue_list.pop())
        d['v_put'] = vectorize(d['queue'].put_nowait)
        self.__dict__.update(d)
        self.__dict__['v_get'] = vectorize(self.get_val)