# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
basedir = os.path.abspath(os.path.join(os.path.dirname('HAOT'), '..'))
sys.path.insert(0, os.path.join(basedir, 'haot', 'src'))
from importlib.metadata import metadata

project = 'haot'
copyright = '2024, Martin E. Liza'
author = 'Martin E. Liza'
release = metadata('haot').get('version')

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',          # Automatically generate API docs
        'sphinx.ext.napoleon',         # Support for NumPy/Google style docstrings
        'sphinx.ext.viewcode',         # Add links to source code
        'sphinx_autodoc_typehints',    # Type hinting in documentation
        'sphinx.ext.mathjax',          # Render LaTeX in HTML output
        ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
autodoc_mock_imports = ['setup']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']

latex_elements = {
    'papersize': 'a4paper',  # Set the paper size to A4 (or 'letterpaper' for US Letter)
    'pointsize': '10pt',     # Set the font size
    'preamble': r'''
\usepackage{amsmath}  % Add math packages
\usepackage{amssymb}  % Additional symbols
''',
}
