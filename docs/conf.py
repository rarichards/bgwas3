# Project information
project = 'bgwas3'
copyright = '2019, Gregory Leeman'
author = 'Gregory Leeman'
release = '0.0.1'

# General configuration
master_doc = 'index'
extensions = [
    'sphinxcontrib.mermaid'
]
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Options for HTML output
html_theme = 'alabaster'
html_static_path = ['_static']
