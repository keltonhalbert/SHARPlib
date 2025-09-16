# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import nwsspc.sharp.calc as sharplib 
import datetime
import breathe
import sys, os 

sys.path.insert(0, os.path.abspath("../"))

version_tuple = sharplib.__version_tuple__
version = sharplib.__version__
year = datetime.datetime.now(datetime.UTC).year
doc_version = os.environ.get('DOC_VERSION', 'dev' if 'dev' in version else version).lstrip('v')

project = 'SHARPlib'
copyright = f'2024-{year}, NOAA/NWS/NCEP Storm Prediction Center'
author = 'Kelton Halbert'
release = f'{version_tuple[0]}.{version_tuple[1]}.{version_tuple[2]}'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.graphviz',
    'sphinx.ext.intersphinx',
    'sphinx.ext.githubpages',
    'breathe',
]

html_static_path = ['_static']
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
add_function_parentheses = False

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

napoleon_numpy_docstring = True 
napoleon_google_docstring = False
napoleon_preprocess_types = True

breathe_projects = {
    "SHARPlib": "doxyxml/"
}
breathe_domain_by_extension = {"h": "cpp"}

breathe_projects_source = {
    "SHARPlib" : (
        "../", [
            "include/SHARPlib/constants.h", 
            "include/SHARPlib/interp.h",
            "include/SHARPlib/layer.h",
            "include/SHARPlib/parcel.h",
            "include/SHARPlib/thermo.h",
            "include/SHARPlib/winds.h",
            "include/SHARPlib/params/convective.h",
            "include/SHARPlib/params/fire.h",
            "include/SHARPlib/params/winter.h",
            "src/SHARPlib/interp.cpp",
            "src/SHARPlib/layer.cpp",
            "src/SHARPlib/parcel.cpp",
            "src/SHARPlib/thermo.cpp",
            "src/SHARPlib/winds.cpp",
            "src/SHARPlib/params/convective.cpp",
            "src/SHARPlib/params/fire.cpp",
            "src/SHARPlib/params/winter.cpp",
        ] 
    )
}

breathe_default_project = "SHARPlib"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "show_toc_level": 2,
    "logo": {
        "text": f"SHARPlib",
        "image_light": "_static/logo.png",
        "image_dark": "_static/logo.png",
    },
    "collapse_navigation": True,
    'external_links': [
        {'name': 'Examples', 'url': 'https://github.com/keltonhalbert/SHARPlib/tree/main/examples/Python'},
        {'name': 'Release Notes', 'url': 'https://github.com/keltonhalbert/SHARPlib/releases'},
    ],
    'icon_links': [
        {
            'name': 'GitHub',
            'url': 'https://github.com/keltonhalbert/SHARPlib',
            'icon': 'fa-brands fa-github',
            'type': 'fontawesome',
        },
        {
            'name': 'Bluesky',
            'url': 'https://bsky.app/profile/stormscale.io',
            'icon': 'fa-brands fa-bluesky',
            'type': 'fontawesome',
        },
    ],
    'use_edit_page_button': False,
    'navbar_align': 'left',
    'navbar_start': ['navbar-logo', 'version-switcher'],
    'navbar_center': ['navbar-nav'],
    'header_links_before_dropdown': 6,
    'navbar_persistent': ['search-button'],
    'navbar_end': ['navbar-icon-links', 'theme-switcher'],
    'switcher': {
        'json_url': 'https://keltonhalbert.github.io/SHARPlib/versions.json',
        'version_match': doc_version
    },
}


html_title = f"{project} v{version} Manual"
html_last_updated_fmt = '%b %d, %Y'
html_context = {"default_mode": "dark"}
html_use_modindex = True
html_copy_source = False
html_domain_indices = False
html_file_suffix = '.html'


intersphinx_mapping = {
    'numpy': ('https://numpy.org/doc/stable/', None),
}
