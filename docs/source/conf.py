import os
import sys 

root_dir = os.path.abspath("../..")   
sys.path.insert(0, root_dir)

project = "TMFC_denoise"
copyright = "2025, Ruslan Masharipov"
author = "Ruslan Masharipov"

release = "v1.4.4"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
root_doc = "index"

extensions = ["sphinx.ext.autodoc", "sphinxcontrib.cairosvgconverter"]

html_theme = 'sphinx_rtd_theme'
html_static_path = ["_static"] 

language = "en"

latex_engine = "pdflatex"   
latex_use_xindy = True   
       
latex_elements = {
    'papersize': 'a4paper',
    'classoptions': ',openany,oneside',
}