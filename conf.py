# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Solar Orbiter Python'
copyright = '2021, David Stansby'
author = 'David Stansby'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.intersphinx',
              'sphinx_gallery.gen_gallery',
              'sphinx_rtd_theme']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

nitpicky = True

# Extension config
# Theme config
# html_theme_options = {'style_nav_header_background': '#ed1b2f'}

# Sphinx gallery
sphinx_gallery_conf = {
    'examples_dirs': 'examples',
    'gallery_dirs': '_examples_build',
    'filename_pattern': r'.*\.py',
    'reference_url': {},
}

# Intersphinx
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'aiapy': ('https://aiapy.readthedocs.io/en/stable', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
    'heliopy': ('https://docs.heliopy.org/en/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable', None),
    'reproject': ('https://reproject.readthedocs.io/en/stable/', None),
    'sunpy': ('https://docs.sunpy.org/en/stable/', None),
    'plasmapy': ('https://docs.plasmapy.org/en/stable', None)
}
