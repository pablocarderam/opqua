# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'Opqua'
copyright = '2023, Pablo Cárdenas'
author = 'Pablo Cárdenas'
release = '1.0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_parser", 
              "sphinx.ext.autodoc",
              "sphinx.ext.intersphinx",
              "sphinx.ext.mathjax",
              "sphinx.ext.viewcode",
              "sphinx.ext.napoleon",
              "sphinx.ext.autosummary",
              # "autoapi.extension",
              "nbsphinx",
              ]

# autoapi_dirs = ['../opqua']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The master toctree document.
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'furo'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
