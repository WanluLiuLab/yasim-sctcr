"""
Configuration file for the Sphinx documentation builder.
"""

# pylint: disable=wrong-import-position, invalid-name

import os

import tomli
from docutils.parsers.null import Parser as NullParser
from sphinx.application import Sphinx

import yasim


def setup(app: Sphinx):
    app.add_source_parser(NullParser)


os.environ['SPHINX_BUILD'] = '1'  # Disable chronolog.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)

# -- Project information -----------------------------------------------------

with open(os.path.join(ROOT_DIR, "pyproject.toml"), "rb") as reader:
    parsed_pyproject = tomli.load(reader)

project = parsed_pyproject["project"]["name"]
author = "&".join([author["name"] for author in parsed_pyproject["project"]["authors"]])
copyright_string = f'2022-2023, {author}'
release = yasim.__version__

# -- General configuration ---------------------------------------------------

html_theme = 'furo'
extensions = [
    # 'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    # 'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    "sphinx.ext.viewcode",
    'myst_nb',
    'sphinx_copybutton',
    'sphinx_design'
]
myst_enable_extensions = ["deflist", "dollarmath"]
exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store',
    '.virtualenv/**'
]

# html_static_path = ['_static']

# Source code suffixes
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'myst-nb',
    '.ipynb': 'null'
}
nb_custom_formats = {
    ".ipynb.py": ["jupytext.reads", {"fmt": "py:percent"}],
    ".ipynb.md": ["jupytext.reads", {"fmt": "md:myst"}]
}

# Insert both docstring of the class and constructor.
autodoc_default_options = {
    'special-members': '__init__',
}

# myst-nb settings
nb_execution_timeout = 1200
nb_execution_mode = "cache"
nb_merge_streams = True
nb_execution_allow_errors = True
