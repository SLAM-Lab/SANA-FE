# Copyright (c) 2026 - The University of Texas at Austin
#  This work was produced under contract #2317831 to National Technology and
#  Engineering Solutions of Sandia, LLC which is under contract
#  No. DE-NA0003525 with the U.S. Department of Energy.
"""
Model registry for the SANA-FE backend.
"""

model_registry = {
    ## Input spike generator ##
    'IN':  {'graph_type': 'input',
            'hw_type': 'soma',
            'sanafe_model': 'input'},

    ## Default Sango Leaky integrate-and-fire neuron ##
    'LIF': {'graph_type': 'neuron',
            'hw_type': 'soma',
            'sanafe_model': 'sango_lif'},

    ## Post-synaptic potential synapse ##
    'PSP': {'graph_type': 'synapse',
            'hw_type': 'synapse',
            'sanafe_model': 'current_based',
            'attribute_types': {'delay': int}},
}
