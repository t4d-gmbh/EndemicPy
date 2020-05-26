from Queue import Empty
from Queue import Queue as SimpleQueue
import sys
import numpy.random as n_rand
from numpy import ndarray
from numpy import array, float64

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


def get_value(kwargs, parameter):
    if parameter in kwargs:
        return kwargs[parameter]
    else:
        return set_default(parameter)


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
    def __init__(
            self,
            distribution_type='exp',
            **kwargs
    ):
        self.distribution_type = distribution_type
        self.reset_random(**kwargs)

    def reset_random(self, **kwargs):
        """
        (Re)set random number generator along with the distribution and queue.

        """
        self._init_draw_fct(**kwargs)
        self.queue = SimpleQueue(maxsize=self._dist_params['size'])
        # fill the queue
        self.fillup()

    def v_put(self, new_values):
        for n_val in new_values:
            self.queue.put_nowait(n_val)

    def v_get(self, shape_obj):
        if isinstance(shape_obj, ndarray):
            length = shape_obj.size
        else:
            length = len(shape_obj)
        return [self.get_val() for _ in xrange(length)]

    def _init_draw_fct(self, **kwargs):
        """
        Initiate the random number generator and the distribution to draw from.
        """
        if self.distribution_type == 'exp':
            self.draw_fct = n_rand.exponential
            self._dist_params = {
                'scale': get_value(kwargs, 'scale'),
                'size': get_value(kwargs, 'size')
                }
            if self._dist_params[
                    'scale'] == 0 or self._dist_params[
                            'scale'] is None:
                self.draw_fct = inf_time
        elif self.distribution_type == 'normal':
            self.draw_fct = n_rand.normal
            self._dist_params = {
                'scale': get_value(kwargs, 'scale'),
                'size': get_value(kwargs, 'size'),
                'loc': get_value(kwargs, 'loc')
                }
            if self._dist_params['scale'] == 0:
                # TOIMPLEMENT
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

    def __getstate__(self):
        d = dict(self.__dict__)
        queue = d.pop('queue')
        simple_queue_list = []
        while True:
            try:
                simple_queue_list.append(queue.get_nowait())
            except Empty:
                break
        d['simple_queue_list'] = simple_queue_list
        # now rebuild the queue for current object
        self._init_draw_fct(**self._dist_params)
        self.queue = SimpleQueue(maxsize=self._dist_params['size'])
        for _el in simple_queue_list:
            self.queue.put_nowait(_el)
        return d

    def __setstate__(self, d):
        if 'simple_queue_list' in d:
            event_queue_list = d.pop('simple_queue_list')[::-1]
            d['queue'] = SimpleQueue(maxsize=d['_dist_params']['size'])
            while len(event_queue_list):
                d['queue'].put_nowait(event_queue_list.pop())
        self.__dict__.update(d)
        self._init_draw_fct(**self._dist_params)
