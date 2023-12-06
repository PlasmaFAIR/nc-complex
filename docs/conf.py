# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import argparse
from breathe import apidoc
import os
import subprocess
import sys

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "nc-complex"
copyright = "2023, Peter Hill"
author = "Peter Hill"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["breathe", "myst_parser"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The default role for text marked up `like this`
default_role = "any"

# Tell sphinx what the primary language being documented is.
primary_domain = "c"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "c"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]

html_context = {
    "github_user": "PlasmaFAIR",
    "github_repo": "nc-complex",
    "github_version": "main",
    "doc_path": "docs",
}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = dict(
    # analytics_id=''  this is configured in rtfd.io
    # canonical_url="",
    repository_url="https://github.com/PlasmaFAIR/nc-complex",
    repository_branch="main",
    path_to_docs="docs",
    use_edit_page_button=True,
    use_repository_button=True,
    use_issues_button=True,
    home_page_in_toc=False,
)

# -- Running doxygen/breathe -------------------------------------------------

breathe_projects = {"nc-complex": "_doxygen/xml"}
breathe_default_project = "nc-complex"
breathe_projects_source = {"nc-complex": ("..", ["include/nc_complex/nc_complex.h"])}
breathe_doxygen_config_options = {
    "PREDEFINED": "NC_COMPLEX_EXPORT=",
    "MACRO_EXPANSION": "YES",
    "ENABLE_PREPROCESSING": "YES",
    "EXPAND_ONLY_PREDEF": "YES",
}
breathe_domain_by_extension = {
    "h": "c",
}
breathe_doxygen_aliases = {
    "rstref{1}": r"\verbatim embed:rst:inline :any:`\1` \endverbatim"
}
