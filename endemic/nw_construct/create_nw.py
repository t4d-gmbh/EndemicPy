__author__ = 'Jonas I Liechti'
#this will be used to create a command line tool, so it is of no importance here.
from argparse import ArgumentParser

from endemic.nw_construct import Graph


def main():
    parser = ArgumentParser(description='Generate a network following a probabilistic framework.',
                                     epilog="Note: Only the options specific for the chosen network_type need to be\
                                      entered."
    )
    parser.add_argument('n', metavar='n', type=int,
                        help='the number of nodes in the network.')
    parser.add_argument('--type', dest='network_type', type=str,
                        help='chose the type of network to construct.')
    parser.add_argument('--method', dest='method', type=str, default='stub',
                        help="set the network construction method (stub or proba, default: stub).")
    parser.add_argument('--file', dest='filename', type=str,
                        help='The file where the network will be stored')
    parser.add_argument('--format', dest='fileformat', type=str, default='edges',
                        help="Specify the format of the destination file (default: edges). If the format is given with \
                        the filename, this argument is ignored.")
    parser.add_argument('--shape', dest='shape', type=float,
                        help='The shape parameter of a distribution (for gamma, exponential and weibull)')
    #is a for weibull
    parser.add_argument('--centre', dest='loc', type=float,
                        help='Centre of the normal distribution.')
    parser.add_argument('--scale', dest='scale', type=float,
                        help='The scale parameter of a distribution (for gamma, normal and exponential)')
    parser.add_argument('--trials', dest='n', type=int,
                        help='The number of trials for a binomial distribution')
    parser.add_argument('--p', dest='p', type=float,
                        help='The success probability (for binomial and geometric).')
    parser.add_argument('--lambda', dest='lam', type=float,
                        help='Expected value for the poisson distribution')
    parser.add_argument('--1+exponent', dest='a', type=float,
                        help='Set the exponent for a power distribution (the dist will be p(x,a)=ax^{a-1}).')
    parser.add_argument('--l', dest='l', type=int,
                        help='In l-partition network, the number of partitions in an l-partition graph.')
    parser.add_argument('--avg_degree', dest='avg_degree', type=float,
                        help="In l-partition network, the expected average degree of each node in the l-partition \
                        network.")
    parser.add_argument('--density_ratio', dest='density_ratio', type=float,
                        help="In l-partition network, the egde-density ratio between inter- and intra-partition links \
                        (eg. 2: inside twice as dense as outside)")
    parser.add_argument('--p_in', dest='p_in', type=float,
                        help='In l-partition network, the connection probability of any pair inside a partition.')
    parser.add_argument('--p_out', dest='p_out', type=float,
                        help='In l-partition network, the connection probability of any pair in-between partitions.')
    #parser.add_argument('integer', metavar='n', type=str, nargs='+',
    #                    help='an integer for the accumulator')
    #parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                    const=sum, default=max,
    #                    help='sum the integers (default: find the max)')

    args = parser.parse_args()
    as_dict = args.__dict__
    filename = as_dict.pop('filename')
    fileformat = as_dict.pop('fileformat')
    if as_dict['network_type'] == 'weibull':
        as_dict['a'] = as_dict.pop('shape')
    if filename is None:
        raise IOError('Please specify a destination to save the graph (--file FILENAME)')
    for key in as_dict.keys():
        if as_dict[key] is None:
            as_dict.pop(key)
    g = Graph(n=as_dict.pop('n'), method=as_dict.pop('method'), **as_dict)
    g.export_graph(filename=filename, fileformat=fileformat)
    return 0
