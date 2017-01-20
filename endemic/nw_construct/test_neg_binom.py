from GraphConstructor import Graph

n=1000
k=4
var=2.0

def get_p(k, var):
    """ 
        In the numpy implementation the probability of success is needed, so
        this returns to probability of success.
    """
    return k/float(var)
def get_r(k, var):
    return (k**2/var)*(1/(var/float(k)-1))
m_g = Graph(
        N=n,
        method='stub',
        p=get_p(k,var),
        n=get_r(k, var)
        )
