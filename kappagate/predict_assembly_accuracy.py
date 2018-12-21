"""This application is experimental."""

import itertools
import networkx
import tatapov
from topkappy import (KappaAgent, KappaSiteState, KappaRule, KappaModel,
                      snapshot_agent_nodes_to_graph)

def linear_graph_to_nodes_list(graph, node_name=None):
    start, end = [n for n in graph.nodes if graph.degree[n] == 1]
    linear_path_nodes = list(networkx.all_simple_paths(graph, start, end))[0]
    if node_name is None:
        return tuple(linear_path_nodes)
    return tuple([graph.nodes[n][node_name] for n in linear_path_nodes])

def predict_assembly_accuracy(slots=None, overhangs=None, duration=1000,
                              initial_quantities=1000, corrective_factor=1.0,
                              annealing_data=('25C', '01h')):
    if overhangs is not None:
        overhangs = list(overhangs)
        slots_overhangs = zip(['LEFT'] + overhangs, overhangs + ['RIGHT'])
        slots = [
            ("p%03d" % i, left, right)
            for i, (left, right) in enumerate(slots_overhangs)
        ]
        
    if isinstance(annealing_data, tuple):
        ad_temp, ad_duration = annealing_data
        annealing_data = tatapov.annealing_data[ad_temp][ad_duration]
    agents = [
        KappaAgent(pos, (left, right))
        for pos, left, right in slots
    ]
    if isinstance(initial_quantities, int):
        initial_quantities = {a: initial_quantities for a in agents}
    all_overhangs = set([o for (pos, l, r) in slots for o in (l, r)
                    if set(o) <= set('ATGC')])
    rules = []
    for agent1, agent2 in itertools.product(agents, agents):
        a1_left, a1_right = agent1.sites
        a2_left, a2_right = agent2.sites
        for site1, site2, a2_side in ((a1_left, a2_left, 'left'),
                                      (a1_left, a2_right, 'right')):
            if (site1 not in all_overhangs) or (site2 not in all_overhangs):
                continue
            if a2_side == 'left':
                ov1, ov2 = site1, site2
            else:
                ov1, ov2 = site1, tatapov.reverse_complement(site2)
            rev_ov1 = tatapov.reverse_complement(ov1)
            rev_ov2 = tatapov.reverse_complement(ov2)
            
            rate = annealing_data[ov1][ov2] + annealing_data[rev_ov1][rev_ov2]
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
    model = KappaModel(
        agents = agents,
        rules = rules,
        initial_quantities=initial_quantities,
        duration=duration,
        snapshot_times={'end': duration}
    )
    simulation_results = model.get_simulation_results()
    expected_slots_order = tuple(pos for pos, _, _ in slots)
    first_slot, last_slot = expected_slots_order[0], expected_slots_order[1]
    end_agents = simulation_results['snapshots']['end']['snapshot_agents']
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