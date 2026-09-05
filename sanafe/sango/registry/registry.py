# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
Baseline model registry for the SANA-FE backend.
"""

model_registry = {
    ## Input spike generator ##
    'IN':  {'graph_type': 'input',
            'hw_type': 'soma',
            'model_type': 'input',
            'param': {},
            'state': {}},

    ## Leaky integrate-and-fire neuron ##
    'LIF': {'graph_type': 'neuron',
            'hw_type': 'soma',
            'model_type': 'sango_lif',
            'param': {},
            'state': {'voltage':   {'dsl': 'voltage',   'default': 0.0},
                      'threshold': {'dsl': 'threshold', 'default': 1.0},
                      'reset':     {'dsl': 'reset',     'default': 0.0},
                      'bias':      {'dsl': 'bias',      'default': 0.0},
                      'leak':      {'dsl': 'leak',      'default': 1.0}}},

    ## Post-synaptic potential synapse ##
    'PSP': {'graph_type': 'synapse',
            'hw_type': 'synapse',
            'model_type': 'current_based',
            'param': {},
            'state': {'delay':    {'dsl': 'delay',  'default': 1.0, 'rep': 'tick'},
                      'weight':   {'dsl': 'weight', 'default': 1.0}}},

    'pLIF': {'graph_type': 'neuron',
             'hw_type': 'soma',
             'model_type': 'sango_plif',
             'hw_name': 'sango_plif',
             'param': {},
             'state': {'voltage':   {'dsl': 'voltage',   'default': 0.0},
                       'threshold': {'dsl': 'threshold', 'default': 1.0},
                       'reset':     {'dsl': 'reset',     'default': 0.0},
                       'bias':      {'dsl': 'bias',      'default': 0.0},
                       'leak':      {'dsl': 'leak',      'default': 1.0},
                       'prob':      {'dsl': 'prob',      'default': 1.0}}},

}
