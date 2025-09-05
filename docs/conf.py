# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import nwsspc.sharp.calc as sharplib 
import breathe
import sys, os 

sys.path.insert(0, os.path.abspath("../"))

project = 'SHARPlib'
copyright = '2025, Kelton Halbert'
author = 'Kelton Halbert'
release = '1.0.3'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.graphviz',
    'breathe',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autoclass_content = "both"
autodoc_member_order = "bysource"
autodoc_docstring_signature = True
autodoc_typehints = "both"
autodo_type_aliases = {
    "PressureLayer": "nwsspc.sharp.calc.layer.PressureLayer",
    "HeightLayer": "nwsspc.sharp.calc.layer.HeightLayer",
    "LayerIndex": "nwsspc.sharp.calc.layer.LayerIndex",
}
autodoc_typehints_format = "short"
autodoc_preserve_defauls = True

napoleon_numpy_docstring = False 
napoleon_google_docstring = True
napoleon_preprocess_types = True

breathe_projects = {
    "SHARPlib": "doxyxml/"
}

breathe_projects_source = {
    "SHARPlib" : (
        "../include/SHARPlib", [
            "constants.h", 
            "interp.h",
            "layer.h",
            "parcel.h",
            "thermo.h",
            "winds.h",
            "params/convective.h",
            "params/fire.h",
            "params/winter.h",
        ] 
    )
}

breathe_default_project = "SHARPlib"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

