project = 'SANA-FE'
author = 'James Boyle'
release = '2.0.20'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
]

html_theme = 'sphinx_rtd_theme'
html_logo = '../logo.svg'

# Add custom CSS
html_static_path = ['_static']
def setup(app):
    app.add_css_file('make_logo_white.css')
