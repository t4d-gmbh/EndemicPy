__author__ = 'Jonas I Liechti'

"""
    This module contains useful functions specific to the model it is in
"""
import os
import sys
import pickle
import glob
import errno
from matplotlib import pyplot as plt
sys.path.append(os.path.abspath('../'))
from nw_construct import Graph
from nw_spread import *

#
# Directories
#

# get the dirs for associated data
resources_path = '../../project_data/'
task = 'gamma_dist'
mode = 'testing'  # leave empty for the real simulations
res_dir = os.path.join(resources_path, task, mode)
pickles_dir = os.path.join(res_dir, 'pickles')
pickles_raw = os.path.join(pickles_dir, 'raw')
pickles_inter = os.path.join(pickles_dir, 'inter')
pickles_res = os.path.join(pickles_dir, 'res')
results_dir = os.path.join(res_dir, 'results')
plots_dir = os.path.join(results_dir, 'plots')

#
# Parameters
#

# pickle names
param_match_pickle = os.path.join(pickles_inter, 'param_match.p')  # contains matching between parameter name and pos.
# Note: pickled results are always named in the *_vs_*.py scripts.
param_comb_pickle = os.path.join(pickles_inter, 'existing_param_comb')  # list of all existing param combinations.

# Stoptime of the selection phase
selection_limit = 400
# Set the names for the strains
name1 = 'wild_type'
wt_id = 0
name2 = 'mutant_1'
mt_id = 1
strains = [name1, name2]

#
# Functions
#
def assure_exist(a_folder):
    """
    Calling this function will check if the folder provided in the argument exists and if not create it.
    :param a_folder: absolute or relative path to the folder to create/check existence
    :type a_folder: str
    :return:
    """
    try:
        os.makedirs(a_folder)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
