__author__ = 'Jonas I Liechti'

DESC = """
    script to run a simple spreading (sis) on an ant colony.
    For this to work you need either to install the endemic package v1.0.0:
        https://tb-git.usys.ethz.ch/j-i-l/EndemicPy/repository/archive.tar.gz?ref=v0.1.0
    Or have the endemic package in the same folder as your script
    (or some subfolder).
"""
import os, sys
from argparse import ArgumentParser
from endemic.nw_construct import TemporalGraph
from endemic.nw_spread import Scenario, ContactSequence, Strain, Pathogen
from random import choice


def job(data_file, beta, mu, starter, output_file, nbr_reps, dt, do):
    """

    Parameter:
    ----------

    :param data_file: File you want to import
    :type data_file: str
    :param beta: Transmission rate
    :type beta: float
    :param mu: Recover rate
    :type mu: float
    :param starter: From where should the infection start. This can either be
        the id of a node or one of the following strings: 'single_random',
        'different_random'
    :type starter: int, float
    :param output_file: The name of the file to where the output should be
        written.
    :type output_file: str
    :param nbr_reps: How many simulations should be run
    :type nbr_rebps: int
    :param dt: Timestep between after which the status of the simulation
        should be logged.
    :param do: optional argument to specify a file in which the detailed output
        of the simulation should be written
    :return:
    """
    print 'Script is running...'
    rec_rate = mu
    # set the recovery rate
    trans_rate = beta
    # and the transmission rate
    # create the network of hosts
    print '- Creating the temporal network...'
    my_graph = TemporalGraph(
        source=data_file,
        start_tag='Starttime',
        stop_tag='Stoptime',
        node1_tag='Tag1',
        node2_tag='Tag2',
        string_values=['Tag1', 'Tag2', 'Box'],
        delimiter=','
    )
    print '\tdone!'
    # use the network to define a contact network
    susceptible_states = 1
    # set all nodes as susceptible in the beginning
    # (see ContactNetwork.__init__ for more info).

    my_contact_sequence = ContactSequence(
            temporal_graph=my_graph, susceptible=susceptible_states
            )
    print my_contact_sequence.all_nodes
    # now that we have the contact network let's define the pathogen
    # We start with creating a pathogen strain.
    wt = Strain(
        name='wild_type', transmission_rate=trans_rate, recover_rate=rec_rate,
    )
    # For this example we simply use a single strain (no competition or other
    # things) so we can go ahead and define the pathogen (which is just a
    # collection of strains really)
    my_pathogen = Pathogen(strains=[wt])

    other_kwargs = {}
    # if dt is provided, pass it along.
    if dt is not None:
        other_kwargs['dt'] = dt
    # we want to have only a maximum number of 1 transmission per contact
    other_kwargs['single_transmission'] = True
    # now that we have the network and the pathogen, we want to put them
    # together and create a spreading scenario:
    my_scenario = Scenario(
        contact_structure=my_contact_sequence,
        pathogen=my_pathogen, **other_kwargs
    )

    # Fro now, no simulation has been done. We were simply initializing our
    # scenario.
    # Now we can start to run simulations.

    # to start let's set stop time for the simulation
    simulation_start = my_graph.t_start
    simulation_end = my_graph.t_stop

    # a scenario can work through phases. Each phase tells the scenario what
    # it should do.
    # So let's create two simple phases our scenario should process:

    # The first phase is just the artificial infection of a single individual:
    if starter == 'random' or starter == 'single_random':
        starter_indivs = 'random'
    elif starter == 'different_random':
        # Here we choose randomly a starting individual for each run and we
        # make sure not to chose an individual twice.
        n = my_scenario.contact_structure.n
        if nbr_reps > n:
            how_many = nbr_reps / n
            starter_indivs = range(n) * how_many
            rest_nodes = range(n)
            for _ in xrange(nbr_reps - how_many * n):
                a_node = choice(rest_nodes)
                rest_nodes.remove(a_node)
                starter_indivs.append(a_node)
        else:
            the_nodes = range(n)
            starter_indivs = []
            for _ in xrange(nbr_reps):
                a_node = choice(the_nodes)
                the_nodes.remove(a_node)
                starter_indivs.append(a_node)
    else:
        # Here we specified the id of the individual which we want to use as an
        # initial seed.
        the_id = my_graph.all_nodes.index(int(starter))
        starter_indivs = [the_id]

    # The second phase is the actual simulation:
    if do:
        spreading_phase = {
                't_start': simulation_start, 't_stop': simulation_end,
                'explicit': True, 'incremental': do
                }
    else:
        spreading_phase = {'t_stop': simulation_end, 'explicit': True}
    # {'t_start': simulation_start, 't_stop': simulation_end}
    # providing t_stop is enough, t_start is set to be my_scenario.t which is 0
    # Here we simply provide the start and the stop time. What the scenario
    # will do is checking who is infected at 't_start'
    # which will be a randomly chose individual (from the initiation phase)
    # and will simulate the spreading of the disease until the scenario reaches
    # the time defined in 't_stop'

    # Now that we have defined the phases we want our scenario to run through,
    # we only need to tell the scenario to run through them:
    print '- Rrunning the simulation...'
    times = []
    counts = []
    for rep in xrange(nbr_reps):
        print '\t %s run.' % str(rep)
        my_scenario.reset()
        if starter == 'different_random':
            print '\tThe initially infected individual is/are: %s' % (
                starter_indivs[rep]
            )
            initiation_phase = {'new_infection': {
                'wild_type': [starter_indivs[rep]]
                }
            }
        else:
            initiation_phase = {
                    'new_infection': {
                        'wild_type': starter_indivs
                        }
                    }
            # This reads: In this phase we want a new infection.
            # The wild_type strain should be used to infect either a random a
            # specific or a set of individuals.

            print '\tThe initially infected individual is/are: %s' % (
                    starter_indivs
                    )
        my_scenario.run([
            initiation_phase,
            spreading_phase
        ])
        print '\t- done!'

        _times = my_scenario.log.keys()
        _times.sort()
        times.append(_times)
        counts.append([my_scenario.log[a_time].count(0) for a_time in _times])

        # At this point, the simulation is done.
        # Now we probably want to access the output.
        # The scenario holds several attributes we can access to get
        # information about its current state and the phases it has been
        # running through.

    print 'All repetitions are done.\n- Writing to output file...'
    with open(output_file, 'w') as f:
        f.write('#time, count\n')
        for j in xrange(len(times)):
            f.write('-- rep %s --\n' % str(j))
            for i in xrange(len(times[j])):
                f.write('%s, %s\n' % (
                    str(times[j][i]), str(counts[j][i]))
                    )


def main():
    parser = ArgumentParser(
        description="""
        This script runs a single SIS simulation on a temporal interaction
        network.

        For this to work you'll need to install the python numpy package and
        the endemic package v0.1.0. How to install NumPy will tell you the
        internet. The endemic package can be found here:

        https://tb-git.usys.ethz.ch/j-i-l/EndemicPy/repository/archive.tar.gz?ref=v0.1.0

        You can either install the package (terminal: python setup.py install)
        or simply copy the endemic package (the folder) in the same folder as
        this script.
        """,
        epilog="PS. Have a nice day!"
    )
    parser.add_argument(
            '-i', '--interactions', dest='data_file', type=str,
            help='Give the path to the interaction data file.'
            )
    parser.add_argument(
            '-b', '--beta', dest='beta', type=float,
            help='This is the per contact transmission rate.'
            )
    parser.add_argument(
            '-m', '--mu', dest='mu', type=float,
            help='This is the recovery rate.'
            )
    parser.add_argument(
            '-s', '--starter', dest='starter', type=str,
            default='single_random', help="""
            Specify how/where the disease should be introduced:

            'single_random': will infect a single random individual,
            'different_random': will choose for each an individual at random
                excluding already chosen individuals.
            """
            )
    parser.add_argument(
            '-o', '--output_file', dest='output_file',
            help='Specify the relative path to the output file.'
            )
    parser.add_argument(
            '-r', '--repetitions', dest='repetitions', type=int, default=1,
            help='Give the number of time the simulation should be repeated.'
            )
    parser.add_argument(
            '-t', '--dt', dest='dt', default=None, type=float,
            help="""
            Specify the time gap with which the current status should be
            reported.
            """
            )
    parser.add_argument(
            '-d', '--detailed_output', dest='do', default='', type=str,
            help="""
            Specify the relative path to the file in which the detailed
            output should be written.
            """
            )
    args = parser.parse_args()
    as_dict = args.__dict__
    data_file = as_dict.pop('data_file')
    beta = as_dict.pop('beta')
    mu = as_dict.pop('mu')
    starter = as_dict.pop('starter')
    output_file = as_dict.pop('output_file')
    reps = as_dict.pop('repetitions')
    dt = as_dict.pop('dt')
    do = as_dict.pop('do')

    job(data_file, beta, mu, starter, output_file, reps, dt, do)


if __name__ == '__main__':
    main()
