# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TMFC_denoise'
copyright = '2025, Ruslan Masharipov'
author = 'Ruslan Masharipov'
release = 'v1.3'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
import os
import sys 

root_dir = os.path.abspath("../..")   
matlab_src_dir = root_dir
sys.path.insert(0, root_dir)


extensions = ["sphinxcontrib.matlab", "sphinx.ext.autodoc"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'

