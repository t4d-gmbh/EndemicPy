<<<<<<< endemic/nw_spread/RateDistribution.py
from numpy import vectorize, array, float64
=======
from numpy import vectorize, random, apply_along_axis
>>>>>>> endemic/nw_spread/RateDistribution.py
from Queue import Empty
from Queue import Queue as SimpleQueue
import sys
from numpy.random import normal, exponential

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


def get_value(kwargs, parameter):
    return kwargs.get(parameter, set_default(parameter))


def docstr_param(*sub):
    def dec(obj):
        obj.__doc.__ = obj.__doc__.format(*sub)
    return dec


@docstr_param(DEFAULT_VALUES['size'], MAX_LIM)
def inf_time(**kwargs):
    """
        Return an array of size {0} with all elements equal to MAX_LIM = {1}.
        Note: You can adapt both the size and the MAX_LIM by passing the
        explicitly as arguments:
            e.g. inf_time(size=10, MAX_LIM=10000)
    """
    return array(
            [MAX_LIM] * kwargs.get('size', DEFAULT_VALUES['size']),
            dtype=float64
            )


@docstr_param(
        DEFAULT_VALUES['size'],
        DEFAULT_VALUES['scale'],
        DEFAULT_VALUES['loc']
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
        :param scale: The scale parameter for the exponential distribution.
            (default=1)
        :param kwargs: Several parameters are possible:
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
    """
    def __init__(self, distribution_type='exp', **kwargs):
        self.distribution_type = distribution_type
        if self.distribution_type == 'exp':
            self.draw_fct = exponential
            self._dist_params = {
                'scale': get_value(kwargs, 'scale'),
                'size': get_value(kwargs, 'size')
                }
        elif self.distribution_type == 'normal':
            self.draw_fct = normal
            self._dist_params = {
                'scale': get_value(kwargs, 'scale'),
                'size': get_value(kwargs, 'size'),
                'log': get_value(kwargs, 'loc')
                }
        else:
            raise self.DistributionTypeError(
                    '\'{0}\' distribution is not implemented yet.'.format(
                        self.distribution_type
                        )
                    )
        # handle the special case of scale == 0
        if not self._dist_params['scale']:
            self.draw_fct = inf_time
        self.queue = SimpleQueue(maxsize=self.size + 1)
        self.v_put = vectorize(self.queue.put_nowait)
        # fill the queue
        self.fillup()
        self.v_get = vectorize(self.get_val)

    class DistributionTypeError(Exception):
        pass

    def fillup(self):
        self.v_put(abs(self.draw_fct(**self._dist_params)))

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
            d['queue'] = SimpleQueue(maxsize=d['size'] + 1)
            while len(event_queue_list):
                d['queue'].put_nowait(event_queue_list.pop())
        self.__dict__.update(d)
        self.__dict__['v_put'] = vectorize(self.queue.put_nowait)
        self.__dict__['v_get'] = vectorize(self.get_val)
        if not self.scale:
            self.queue = SimpleQueue(maxsize=self.size + 1)
            # this is specific to the queue, thus "reinit" here
            self.v_put = vectorize(self.queue.put_nowait)
            self.draw_fct = inf_time

            self.fillup()
