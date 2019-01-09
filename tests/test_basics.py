import os
import matplotlib
matplotlib.use("Agg")
from kappagate import (overhangs_list_to_slots, predict_assembly_accuracy,
                       plot_colony_picking_graph, success_rate_facts,
                       plot_circular_interactions)

def test_basic_success_rate_prediction():
    overhangs=  ['GGAG', 'GGCA', 'TCGC', 'CAGT', 'TCCA',
                 'GAAT', 'AGTA', 'TCTT', 'CAAA', 'GCAC',
                 'AACG', 'GTCT', 'CCAT']
    slots = overhangs_list_to_slots(overhangs)
    success_rate, _, _ = predict_assembly_accuracy(
        slots, initial_quantities=5000)
    assert success_rate > 0.95


    overhangs = ['GGAG', 'GATA', 'GGCA', 'GGTC', 'TCGC',
                 'GAGG', 'CAGT', 'GTAA', 'TCCA', 'CACA',
                 'GAAT', 'ATAG', 'AGTA', 'ATCA', 'TCTT',
                 'AGGT', 'CAAA', 'AAGC', 'GCAC', 'CAAC',
                 'AACG', 'CGAA', 'GTCT', 'TCAG', 'CCAT']
    slots = overhangs_list_to_slots(overhangs)
    success_rate, _, _ = predict_assembly_accuracy(
        slots, initial_quantities=5000)
    assert 0.8 < success_rate < 0.9


    overhangs=  ['GGAG', 'GGTC', 'AGCA', 'CAGT', 'GGTA',
                 'GAAT', 'GGTT', 'TCTT', 'GGTG', 'GCAC',
                 'AGCG', 'GTCT', 'CCAT']
    slots = overhangs_list_to_slots(overhangs)
    success_rate, _, _ = predict_assembly_accuracy(
        slots, initial_quantities=5000)
    assert 0.3 < success_rate < 0.5

def test_plot_circular_interactions():
    overhangs = ['TAGG', 'GACT', 'GGAC', 'CAGC',
                 'GGTC', 'GCGT', 'TGCT', 'GGTA',
                 'CGTC', 'CTAC', 'GCAA', 'CCCT']
    slots = overhangs_list_to_slots(overhangs)
    plot_circular_interactions(
        slots, annealing_data=('25C', '01h'), rate_limit=200)

def test_success_rate_facts():
    overhangs = ['TAGG', 'GACT', 'GGAC', 'CAGC',
                 'GGTC', 'GCGT', 'TGCT', 'GGTA',
                 'CGTC', 'CTAC', 'GCAA', 'CCCT']
    slots = overhangs_list_to_slots(overhangs)
    predicted_rate, _, _ = predict_assembly_accuracy(slots)
    plot_colony_picking_graph(success_rate=predicted_rate)
