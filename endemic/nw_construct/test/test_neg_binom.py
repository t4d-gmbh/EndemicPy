from ..GraphConstructor import Graph

n=100000
k=4
var=2.0
#distro = 'negative_binomial'
distro = 'gamma'

def get_p(k, var):
    """ 
        In the numpy implementation the probability of success is needed, so
        this returns to probability of success.
    """
    return k/float(var)

def get_r(k, var):
    #return (k**2/var)*(1/(var/float(k)-1))
    return int(k**2/(var - k))

def get_scale(k,var):
    return var / float(k)

def get_shape(k,var):
    return k ** 2 / float(var)

distribution_params = {}
distribution_params['network_type'] = distro
if var== 0.0:
    distribution_params['network_type'] = 'uniform'
    distribution_params['degree'] = k
elif distro == 'negative_binomial':
    distribution_params['p'] = get_p(k,var)
    distribution_params['n'] = get_r(k,var)
elif distro == 'gamma':
    distribution_params['scale'] = get_scale(k,var)
    distribution_params['shape'] = get_shape(k,var)

m_g = Graph(
        N=n,
        method='stub',
        **distribution_params
        )
