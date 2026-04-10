# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RNAquarium'
copyright = '2026, Biohub'
author = (
    'Yttria Aniseia, Eric Waltari, Max Frank, Leandro Lima, '
    'Gibraan Rahman, Andy Zhou, Yang-Joon Kim, Hejin Huang, '
    'Yasin Şenbabaoğlu, Duo Peng, and Keir Balla'
)
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- MyST-Parser configuration -----------------------------------------------

myst_enable_extensions = [
    "attrs_block",
    "colon_fence",
]
myst_number_code_blocks = ["bash", "groovy"]
myst_heading_anchors = 3  # auto-generate anchors for h1-h3

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Logo in the upper-left corner of the sidebar
# Place your logo image as docs/_static/logo.png
html_logo = '_static/logo.png'

html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'includehidden': True,
    'titles_only': False,
    'logo_only': False,       # False = show project name alongside logo
}

# Title shown in the browser tab
html_title = 'RNAquarium Documentation'
html_short_title = 'RNAquarium'
