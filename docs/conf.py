import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../bgwas3"))
import bgwas3

project = "bgwas3"
copyright = "2019, Gregory Leeman"
author = "Gregory Leeman"
release = "0.0.1"
master_doc = "index"
extensions = [
    "sphinx.ext.graphviz",
    "sphinx.ext.autodoc",
    "sphinx.ext.graphviz"
]
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "alabaster"
html_static_path = ["_static"]
autodoc_member_order = "bysource"
