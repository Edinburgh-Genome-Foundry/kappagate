"""This application is experimental."""

import itertools
import networkx
import tatapov
from topkappy import (KappaAgent, KappaSiteState, KappaRule, KappaModel,
                      snapshot_agent_nodes_to_graph)

from .tools import overhangs_list_to_slots, linear_graph_to_nodes_list

def slots_to_agents_and_rules(slots, annealing_data=('25C', '01h'),
                              corrective_factor=1.0):
    """Generate Topkappy rules and agents objects modeling parts interactions.
    
    Parameters
    ----------

    slots
      A list [(slot_name, left_overhang, right_overhang), ...]
    
    annealing_data
      Either a pandas dataframe or a couple (temperature, duration) indicating
      an experimental dataset from Potapov et al. 2018
    
    corrective_factor
      A factor that can be applied to decrease (when <1) or increase (>1)
      the differences in affinity in the dataset.
    
    Returns
    -------

    agents, rules
      Lists of Topkappy agents and rules, ready to be fed to a KappaModel
    """
    if isinstance(annealing_data, tuple):
        ad_temp, ad_duration = annealing_data
        annealing_data = tatapov.annealing_data[ad_temp][ad_duration]
    agents = [
        KappaAgent(pos, (left, right))
        for pos, left, right in slots
    ]
    all_overhangs = set([o for (pos, l, r) in slots for o in (l, r)
                    if set(o) <= set('ATGC')])
    rules = []
    for agent1, agent2 in itertools.product(agents, agents):
        _, a1_right = agent1.sites
        a2_left, a2_right = agent2.sites
        for site1, site2, a2_side in ((a1_right, a2_left, 'left'),
                                      (a1_right, a2_right, 'right')):
            if (site1 not in all_overhangs) or (site2 not in all_overhangs):
                continue
            if a2_side == 'left':
                ov1, ov2 = site1, tatapov.reverse_complement(site2)
            else:
                ov1, ov2 = site1, site2
            # rev_ov1 = tatapov.reverse_complement(ov1)
            # rev_ov2 = tatapov.reverse_complement(ov2)
            
            rate = annealing_data[ov1][ov2] #+ annealing_data[rev_ov2][rev_ov1]
            if rate == 0:
                continue
            rules.append(KappaRule(
                '%s-left.%s-%s' % (agent1.name, agent2.name, a2_side),
                [
                    KappaSiteState(agent1.name, site1, '.'),
                    KappaSiteState(agent2.name, site2, '.')
                ],
                '->',
                [
                    KappaSiteState(agent1.name, site1, '1'),
                    KappaSiteState(agent2.name, site2, '1')
                ],
                rate=rate ** corrective_factor
            ))
    return agents, rules

def predict_assembly_accuracy(slots, duration=1000, initial_quantities=1000,
                              corrective_factor=1.0,
                              annealing_data=('25C', '01h')):
    """Predict the accuracy of the assembly (proportion of good clones).
    
    Parameters
    ----------

    slots
      A list [(slot_name, left_overhang, right_overhang), ...]
    
    annealing_data
      Either a pandas dataframe or a couple (temperature, duration) indicating
      an experimental dataset from Potapov et al. 2018
    
    duration
      Virtual duration of the Kappa complexation simulation experiments.
      A large number ensures that the experiment comes to equilibrium, and
      is not necessarily longer. So keep it large
    
    initial_quantities
      Either a dict {slot_name: initial_quantity} or an integer in case all
      agents start the simulation with the same initial quantity. The higher
      the initial quantities, the less noisy the results of the simulation.
    
    corrective_factor
      A factor that can be applied to decrease (when <1) or increase (>1)
      the differences in affinity in the dataset.
    
    Returns
    -------

    proportion, other_constructs, simulation_results
      Where proportion is the proportion of good clones, other_constructs is
      a dict {parts_tuple: proportion} showing the proportion of circular
      constructs (in bad clones), and simulation_results is the topkappy
      simulation results object.
    """
    agents, rules = slots_to_agents_and_rules(
        slots, annealing_data=annealing_data,
        corrective_factor=corrective_factor)
    if isinstance(initial_quantities, int):
        initial_quantities = {a: initial_quantities for a in agents}
    model = KappaModel(
        agents=agents,
        rules=rules,
        initial_quantities=initial_quantities,
        duration=duration,
        snapshot_times={'end': duration}
    )
    simulation_results = model.get_simulation_results()
    expected_slots_order = tuple(pos for pos, _, _ in slots)
    first_slot, last_slot = expected_slots_order[0], expected_slots_order[1]
    snapshots = simulation_results['snapshots']
    end_time = 'end' if 'end' in snapshots else 'deadlock' 
    end_agents = snapshots[end_time]['snapshot_agents']
    filtered_agents = [
        (freq, snapshot_agent_nodes_to_graph(nodes, with_ports=False))
        for freq, nodes in end_agents
        if any(node['node_type'] == first_slot for node in nodes)
        and any(node['node_type'] == last_slot for node in nodes)
    ]
    n_filtered_agents = sum(fa[0] for fa in filtered_agents)
    filtered_agents_with_slots = {
        linear_graph_to_nodes_list(graph, node_name='node_name'):
        1.0 * freq / n_filtered_agents
        for freq, graph in filtered_agents
    }
    score = (filtered_agents_with_slots.get(expected_slots_order, 0) +
             filtered_agents_with_slots.get(expected_slots_order[::-1], 0))
    return score, filtered_agents_with_slots, simulation_results