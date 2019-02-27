import os
import networkx as nx
from Bio import SeqIO
from dnacauldron import (RestrictionLigationMix, autoselect_enzyme,
                         get_overhangs_from_record)
from dnacauldron.tools import reverse_complement
from snapgene_reader import snapgene_file_to_seqrecord

def parts_records_to_slots(parts_records, enzyme='auto'):
    """Return slots from parts records, ready to feed to other methods.
    
    Parameters
    ----------
    parts_records
      A list of Biopython records of the parts of the assembly. Do NOT
      include the backbone.
    
    enzyme
      Name of a Type-2S enzyme or "auto" to select automatically based
      on restriction sites in the records sequences.


    Returns
    -------
    slots
      A list [(slot_name, left_overhang, right_overhang), ...] ready to be fed
      to the other Kappagate methods.

    """


    if enzyme == 'auto':
        enzyme = autoselect_enzyme(parts_records)
    mix = RestrictionLigationMix(parts_records, enzyme=enzyme)
    slots_parts = mix.compute_slots()
    graph = mix.slots_graph(with_overhangs=False)
    slots = [
        (list(slots_parts[(s1, s2)])[0], s1, s2)
        for s1, s2 in linear_graph_to_nodes_list(graph)
    ]
    for i in list(range(len(slots) - 1)):
        name, left, right = slots[i + 1]
        if slots[i][2] != left:
            if slots[i][2] == reverse_complement(right):
                slots[i + 1] = (name, reverse_complement(right),
                                reverse_complement(left))
            else:
                name, left, right = slots[i]
                slots[i] = (name, reverse_complement(right),
                            reverse_complement(left))
    return ([('backbone-left', 'LEFT', slots[0][1])] +
            slots +
            [('backbone-right', slots[-1][2], 'RIGHT')])

def _find_backbone_center(record, backbone_annotations=()):
    """Find an annotation from the backbone, return the index of its center"""
    record.features = [f for f in record.features
                       if f.location is not None]
    for feature in record.features:
        for qualifier in feature.qualifiers.values():
            qualifiers = str(qualifier)
            if any([ann in qualifiers for ann in backbone_annotations]):
                return int((feature.location.start + feature.location.end)/2)
    raise ValueError("Could not find any of the following in record %s: %s" % (
                     record.id, ", ".join(backbone_annotations)))

def construct_record_to_slots(record, backbone_annotations=()):
    """Return slots from a construct record, ready to feed to other methods.
    
    Parameters
    ----------
    record
      A biopython record of an assembly construct, either created by
      DnaCauldron, or with explicit annotations with  feature type "homology"
      to indicate overhangs.
    
    backone_annotations
      Texts that can be found in the annotations located in the "backbone part"
      of the provided record. e.g. ['AmpR', 'Origin'] etc.

    Returns
    -------
    slots
      A list [(slot_name, left_overhang, right_overhang), ...] ready to be fed
      to the other Kappagate methods.

    """
    backbone_center = _find_backbone_center(
        record, backbone_annotations=backbone_annotations)
    overhangs = get_overhangs_from_record(record, with_locations=True)
    if overhangs is None:
        raise ValueError("Could not find any overhang in the provided record "
                         "with id %s" % record.id)
    overhangs = ([o for loc, o in overhangs if loc > backbone_center] +
                 [o for loc, o in overhangs if loc <= backbone_center])
    return overhangs_list_to_slots(overhangs)

def linear_graph_to_nodes_list(graph, node_name=None):
    """Return a list of node names as they appear in the linear graph."""
    start, end = [n for n in graph.nodes if graph.degree[n] == 1]
    linear_path_nodes = list(nx.all_simple_paths(graph, start, end))[0]
    if node_name is None:
        return tuple(linear_path_nodes)
    return tuple([graph.nodes[n][node_name] for n in linear_path_nodes])

def overhangs_list_to_slots(overhangs):
    """Return slots from a list of overhangs, ready to feed to other methods.
    
    Parameters
    ----------
    overhangs
      A list of the form ['ATGC', 'TTGC', 'TTAC'...] of the overhangs as they
      appear in the modelled assembly (i.e. the first part has overhangs
      ATGC-TTGC, the second TTGC-TTAC, etc.)
    
    backone_annotations
      Texts that can be found in the annotations located in the "backbone part"
      of the provided record. e.g. ['AmpR', 'Origin'] etc.

    Returns
    -------
    slots
      A list [(slot_name, left_overhang, right_overhang), ...] ready to be fed
      to the other Kappagate methods.

    """
    overhangs = list(overhangs)
    slots_overhangs = zip(['LEFT'] + overhangs, overhangs + ['RIGHT'])
    overhangs = [
        ("p%03d" % i, left, right)
        for i, (left, right) in enumerate(slots_overhangs)
    ]
    overhangs[0] = ('backbone-left', *overhangs[0][1:])
    overhangs[-1] = ('backbone-right', *overhangs[-1][1:])
    return overhangs
    

def load_record(filename, linear=True, id='auto', upperize=True):
    if hasattr(filename, 'read'):
        record = SeqIO.read(filename, "genbank")
        if id == 'auto':
            raise ValueError("Can't have id == 'auto' when reading filelikes.")
    elif filename.lower().endswith(("gb", "gbk")):
        record = SeqIO.read(filename, "genbank")
    elif filename.lower().endswith(('fa', 'fasta')):
        record = SeqIO.read(filename, "fasta")
    elif filename.lower().endswith('.dna'):
        record = snapgene_file_to_seqrecord(filename)
    else:
        raise ValueError('Unknown format for file: %s' % filename)
    if upperize:
        record = record.upper()
    record.linear = linear
    if id == 'auto':
        id = record.id
        if id in [None, '', "<unknown id>", '.', ' ']:
            id = os.path.splitext(os.path.basename(filename))[0]
            record.name = id.replace(" ", "_")[:20]
        record.id = id
    elif id is not None:
        record.id = id
        record.name = id.replace(" ", "_")[:20]
    return record