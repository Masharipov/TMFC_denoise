import os
import sys 

root_dir = os.path.abspath("../..")   
matlab_src_dir = root_dir
sys.path.insert(0, root_dir)

project = "TMFC_denoise"
copyright = "2025, Ruslan Masharipov"
author = "Ruslan Masharipov"

release = "v1.3"

extensions = ["sphinxcontrib.matlab", "sphinx.ext.autodoc"]
templates_path = ["_templates"]

html_theme = 'sphinx_rtd_theme'

