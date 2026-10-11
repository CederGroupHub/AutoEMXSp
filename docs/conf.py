# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import re
import sys
sys.path.insert(0, os.path.abspath(".."))



# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AutoEMX'
copyright = '2026, Andrea Giunto'
author = 'Andrea Giunto'
# Version of the package, read from autoemx/_version.py (its single source) so that it follows each release
with open(os.path.join(os.path.dirname(__file__), "..", "autoemx", "_version.py"), encoding="utf-8") as _f:
    release = re.search(r'^__version__ = "([^"]+)"', _f.read(), re.M).group(1)
version = release



# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_js_files = ['expand_user_docs.js']
html_logo = '_static/logo/autoemx-logo-dark.svg'
html_favicon = '_static/logo/favicon.ico'

# Render the full navigation tree on every page (with expand/collapse buttons),
# so expand_user_docs.js can open the User Documentation section by default
html_theme_options = {
    'collapse_navigation': False,
    # Logo already spells the project name; navy matches the logo and the GUI header
    'logo_only': True,
    'style_nav_header_background': '#13263d',
}



# -- Options for autodoc -------------------------------------------------
# Show members in the order they appear in source
autodoc_member_order = "bysource"

# Respect __all__ for user docs
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}