import numpy as np
import treecorr
import healpy as hp
import matplotlib.pyplot as plt
import read_functions
import astropy.io.fits as fits
import os

###  MODIFY ONCE WE HAVE A PACKAGE SETUP
import sys, os
current_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(current_dir)

sys.path.append(parent_dir)


from io.use_config import load_yaml_config

load_yaml_config(parent_dir+'/config/t19lite.yml')
