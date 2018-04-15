from numpy import vectorize, array, float64, random, apply_along_axis
from Queue import Empty
from Queue import Queue as SimpleQueue
import sys
import numpy.random as n_rand
from numpy import ndarray

# this is the time it takes if a rate of 0 is given
MAX_LIM = 100000
# this is the default length of the queue holding drawn event times
DEFAULT_VALUES = {
    'scale': 1.0,
    'loc': 0.0,
    'size': 1000
    }


def set_default(parameter):
    sys.stdout.write(
        'INFO: Parameter \'{0}\' was not provided. It is set to its default \
                value of {1}'.format(parameter, DEFAULT_VALUES[parameter])
                )
    sys.stdout.flush()
    return DEFAULT_VALUES[parameter]


def no_mut(dummy, length):
    return [MAX_LIM] * length
    # TODO: alternative
    # return array([MAX_LIM] * length, dtype=float64)


def get_value(kwargs, parameter):
    if parameter in kwargs:
        return kwargs[parameter]
    else:
        return set_default(parameter)
    return kwargs.get(parameter, set_default(parameter))


def inf_time(**kwargs):
    """
        Return an array of size {0} with all elements equal to MAX_LIM = {1}.
        Note: You can adapt both the size and the MAX_LIM by passing the
        explicitly as arguments:
            e.g. inf_time(size=10, MAX_LIM=10000)
    """.format(DEFAULT_VALUES['size'], MAX_LIM)
    return array(
            [MAX_LIM] * kwargs.get('size', DEFAULT_VALUES['size']),
            dtype=float64
            )


class Distro(object):

    """
        This class holds a queue of random times drawn from a given
        distribution with a specified scale.

        Parameters:
        -----------
        :param distribution_type: Specifies which distribution to use.
            Possible values are: 'exp' for exponential and 'norm' for gaussian.
            The default choice is 'exp'.
        :type distribution_type: str
        # :param seed: provide a seed for the numpy random number generator
        :param kwargs: Several parameters are possible:
            - scale: The scale parameter for the exponential distribution.
                (default=1)
            - size: Predefined size of the queue.
                (default: {0})
            - set of parameters specific to the distribution that is
                chosen. For details about the requested parameters check the
                distribution definitions:
                    - exp: numpy.random.exponential
                    - norm: numpy.random.normal
                Note:
                    - If the parameter 'scale' is not given, the default value
                        of scale={1} is chosen.
                    - If the parameter 'loc' is not given, a default value of
                        loc = {2} is chosen.
    """.format(
        DEFAULT_VALUES['size'],
        DEFAULT_VALUES['scale'],
        DEFAULT_VALUES['loc']
        )
    def __init__(self,
            distribution_type='exp',
            # size=DEFAULT_VALUES['size'],
            # seed=None,
            **kwargs):
        # self.size = size
        # self._init_seed(seed)
        self.distribution_type = distribution_type
        self._init_draw_fct(**kwargs)
        # handle the special case of scale == 0
        self.queue = SimpleQueue(maxsize=self._dist_params['size'])
        #self.v_put = vectorize(self.queue.put_nowait)
        # fill the queue
        self.fillup()
        #self.v_get = vectorize(self.get_val)

    def v_put(self, new_values):
        for n_val in new_values:
            self.queue.put_nowait(n_val)

    def v_get(self, shape_obj):
        if isinstance(shape_obj, ndarray):
            length = shape_obj.size
        else:
            length = len(shape_obj)
        return [self.get_val() for _ in xrange(length)]


    # def _init_seed(self, seed):
    #     if seed is None:
    #         n_rand.seed()
    #         self.seed = n_rand.get_state()
    #     else:
    #         assert isinstance(seed, tuple), \
    #                 'seed, if provided, must be a tuple'
    #         n_rand.set_state(seed)
    #         self.seed = n_rand.get_state()

    def _init_draw_fct(self, **kwargs):
        if self.distribution_type == 'exp':
            self.draw_fct = n_rand.exponential
            self._dist_params = {
                'scale': get_value(kwargs, 'scale'),
                'size': get_value(kwargs, 'size')
                }
            if self._dist_params['scale'] == 0:
                self.draw_fct = inf_time
        elif self.distribution_type == 'normal':
            self.draw_fct = n_rand.normal
            self._dist_params = {
                'scale': get_value(kwargs, 'scale'),
                'size': get_value(kwargs, 'size'),
                'loc': get_value(kwargs, 'loc')
                }
            if self._dist_params['scale'] == 0:
                #TOIMPLEMENT
                # We are at a dirak delta
                pass
        else:
            raise self.DistributionTypeError(
                    '\'{0}\' distribution is not implemented yet.'.format(
                        self.distribution_type
                        )
                    )

    class DistributionTypeError(Exception):
        pass

    def fillup(self):
        self.v_put(self.draw_fct(**self._dist_params))
        #self.v_put(abs(self.draw_fct(**self._dist_params)))
        return 0

    def get_val(self, a=None):
        """
        Function returning a value drawn form the associated distribution.
        """
        try:
            return self.queue.get_nowait()
        except Empty:
            self.fillup()
            return self.queue.get_nowait()

    def get_times(limit):
        """
        Return an array of events happening before the limit.
        
        Parameters:
        -----------
            :param limit: upper limit for the event times sequence. Only events
                happening before the limit will be returned.
            :type limit: float, int
        """
        # draw form the dist until we reach the limit time.    
        pass
    
        
    #def v_get(self, an_array):
    #    #return map(self.get_val, xrange(an_array.size))
    #    return apply_along_axis(self.get_val, 0, an_array)
    # to transform the priority queue holding the upcoming events into a pickle-abel list
    

    def __getstate__(self):
        d = dict(self.__dict__)
        queue = d.pop('queue')
        # v_put is not pickle-able
        #del d['v_put']
        #del d['v_get']
        simple_queue_list = []
        while True:
            try:
                simple_queue_list.append(queue.get_nowait())
            except Empty:
                break
        d['simple_queue_list'] = simple_queue_list
        # now rebuild the queue for current object
        # self._init_seed(self.seed)
        self._init_draw_fct(**self._dist_params)
        self.queue = SimpleQueue(maxsize=self._dist_params['size'])
        for _el in simple_queue_list:
            self.queue.put_nowait(_el)
        #self.v_put = vectorize(self.queue.put_nowait)
        #self.v_get = vectorize(self.get_val)
        return d

    # to load the pickled event list back into a priority queue
    def __setstate__(self, d):
        if 'simple_queue_list' in d:
            event_queue_list = d.pop('simple_queue_list')[::-1]
            d['queue'] = SimpleQueue(maxsize=d['_dist_params']['size'])
            while len(event_queue_list):
                d['queue'].put_nowait(event_queue_list.pop())
        self.__dict__.update(d)
        # not sure if the seed for numpy.random is initiated only here
        # self._init_seed(d.get('seed', None))
        self._init_draw_fct(**self._dist_params)
