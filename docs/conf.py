project = 'SANA-FE'
author = 'James Boyle'
release = '2.0.20'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',  # For automatic summary tables
    'sphinx.ext.napoleon',  # For Google style docstrings
]

html_theme = 'sphinx_rtd_theme'
html_logo = '../sana_fe_logo.svg'

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Generate autosummary automatically
autosummary_generate = True

# Add custom CSS
html_static_path = ['_static']
def setup(app):
    app.add_css_file('make_logo_white.css')

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True