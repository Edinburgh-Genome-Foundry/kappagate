import os
import matplotlib
matplotlib.use("Agg")
from kappagate import (overhangs_list_to_slots, predict_assembly_accuracy,
                       plot_colony_picking_graph, success_rate_facts,
                       plot_circular_interactions, load_record,
                       parts_records_to_slots, construct_record_to_slots)
import flametree

records_dict = {
    name: load_record(os.path.join('tests', 'data', 'records', name + '.gb'),
                      id=name, topology='circular')
    for name in ("partA", "partB", "partC", "assembled_construct")
}

def test_basic_success_rate_prediction():
    overhangs =  ['GGAG', 'GGCA', 'TCGC', 'CAGT', 'TCCA',
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
    assert 0.8 < success_rate < 0.92


    overhangs =  ['GGAG', 'GGTC', 'AGCA', 'CAGT', 'GGTA',
                  'GAAT', 'GGTT', 'TCTT', 'GGTG', 'GCAC',
                  'AGCG', 'GTCT', 'CCAT']
    slots = overhangs_list_to_slots(overhangs)
    success_rate, _, _ = predict_assembly_accuracy(
        slots, initial_quantities=5000)
    assert 0.2 < success_rate < 0.4

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
    predicted_rate, _, _ = predict_assembly_accuracy(slots, duration=10)
    plot_colony_picking_graph(success_rate=predicted_rate)

def test_parts_records_to_slots():
    records = [records_dict[n] for n in ["partA", "partB", "partC"]]
    slots = parts_records_to_slots(records, enzyme='auto')
    assert slots == [('backbone-left', 'LEFT', 'ATTG'),
                     ('partA', 'ATTG', 'GGCT'),
                     ('partB', 'GGCT', 'GGGC'),
                     ('partC', 'GGGC', 'GGCA'),
                     ('backbone-right', 'GGCA', 'RIGHT')]
    

def test_construct_record_to_slots():
    record = records_dict['assembled_construct']
    slots = construct_record_to_slots(record, backbone_annotations='receptor')
    assert slots == [('backbone-left', 'LEFT', 'ATTG'),
                     ('p001', 'ATTG', 'GGCT'),
                     ('p002', 'GGCT', 'GGGC'),
                     ('p003', 'GGGC', 'GGCA'),
                     ('backbone-right', 'GGCA', 'RIGHT')]