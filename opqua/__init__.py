"""
Opqua.

An epidemiological modeling framework for population genetics and evolution.
"""

__version__ = "0.2.1"
__author__ = 'Pablo Cardenas'
__credits__ = 'Massachusetts Institute of Technology'

# init file to make this directory into a python package
from os.path import dirname, basename, isfile
import glob
modules = glob.glob(dirname(__file__)+"/*.py")
__all__ = [ basename(f)[:-3] for f in modules if isfile(f)]
