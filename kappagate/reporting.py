import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import itertools

from .predict_assembly_accuracy import slots_to_agents_and_rules
from .tools import overhangs_list_to_slots, linear_graph_to_nodes_list


def plot_colony_picking_graph(success_rate=None, ax=None):
    """Plot a generic graph of colonies to pick in different scenarios.

    Parameters
    ----------

    success_rate
      If provided (betweem 0 and 1), it will be highlighted on the map
    
    ax
      A matplotlib ax (one is created and returned if none is provided)
    
    Returns
    -------

    ax
      The matplotlib ax of the plot
    """
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(15, 5))
    ax.set_title("Number of colonies to pick to get at least one "
                 "good clone with certainty X%\n",
                 fontdict=dict(size=16))
    ax.set_xlabel("Proportion of good clones", fontdict=dict(size=16))
    ax.set_ylabel("Certainty level", fontdict=dict(size=16))
    ax.set_yscale('log')
    # ax.set_xscale('log')

    # ax.invert_yaxis()

    X = np.linspace(0.05, 0.999, 1000)
    Y = np.linspace(0.001, 0.500, 1000)[::-1]
    XX, YY = np.meshgrid(X, Y)
    ZZ = np.log(YY)/np.log(1-XX)
    levels = (list(range(1, 10)) +
              list(range(10, 20, 2)) +
              list(range(20, 30, 5)) +
              list(range(30, 55, 10))) 
    # ax.contourf(XX, YY, ZZ, levels=levels)

    xticks = np.arange(0.1, 1.01, 0.1)
    yticks = [0.5, 0.2, 0.1, 0.05, 0.01]
    ax.set_xticks(xticks)
    ax.set_xticklabels([("%.02f" % (x)).rstrip('0') for x in xticks])
    ax.set_yticks(yticks)
    ax.set_yticklabels(["%d%%" % (100 * (1 - y)) for y in yticks])
    cs = ax.contour(XX, YY, ZZ, levels=levels, colors=['black'])
    ax.clabel(cs, inline=1, fontsize=10, fmt='%d',
              manual=[(1.0 - (0.025) ** (1.0/x), 0.025) for x in levels])
    if success_rate is not None:
        ax.axvline(x=success_rate, ls=':', color='red')
        ax.plot(len(yticks) * [success_rate],
                [y for y in yticks],
                lw=0, marker='o', c='r', markeredgecolor='r',
                markerfacecolor='white')
    return ax

def min_trials_for_one_success(success_rate, certainty):
    """Return the minimal number of trials to be X% certain to have at least
    one success."""
    if success_rate == 1:
      return 1
    if success_rate == 0:
      return np.Inf
    return np.ceil(np.log(1 - certainty) / np.log(1 - success_rate))

def average_trials_until_success(success_rate):
    """Return the average number of trials before a success is encountered."""
    if success_rate == 0:
      return np.Inf
    return 1.0 / success_rate

def success_rate_facts(success_rate, plain_text=True):
    """Return relevant stats for the given success rate.
    
    
    Returns
    -------
    plain_text (if plain_text=True)
      A plain text as follows: "The valid colony rate is 47.7%. Expect 1.9
      clones in average until success. Pick 5 clones or more for 95% chances
      of at least one success. If no success after 8 clones, there is
      likely another problem (p-value=0.01)"

    dict (if plain_text=False)
      Dict containing the same infos as above: dict(success_rate_percent,
      average_colonies, min_trials_q95, max_trials_q99)
        )

    """
    if not plain_text:
        return dict(
            success_rate_percent=100 * success_rate,
            average_colonies=average_trials_until_success(success_rate),
            min_trials_q95=min_trials_for_one_success(success_rate, 0.95),
            max_trials_q99=min_trials_for_one_success(success_rate, 0.99)
        )
    results =  success_rate_facts(success_rate, plain_text=False)
    return (
        "The valid colony rate is %(success_rate_percent).1f%%. Expect "
        "%(average_colonies).1f clones in average until "
        "success. Pick %(min_trials_q95)d clones or more for 95%% chances of "
        "at least one success. If no success after %(max_trials_q99)d clones, "
        "there is likely another problem (p-value=0.01)") % results

def plot_circular_interactions(slots, annealing_data=('25C', '01h'),
                               corrective_factor=1.0, rate_limit=200, ax=None):
    """Plot the slots circularly, show the strength of overhangs interactions.
    
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
    
    rate_limit
      Any interaction with a "rate" (value in the annealing_data) below this
      value with not be shown.
    
    ax
      A matplotlib ax (one is created and returned if none is provided)
    
    Returns
    -------
    ax
      The matplotlib ax of the plot

    """ 
    agents, rules = slots_to_agents_and_rules(
        slots, annealing_data=annealing_data,
        corrective_factor=corrective_factor)


    graph = nx.Graph([((s[0], s[1]), (s[0], s[2])) for s in slots])
    expected_interactions = []
    for n1, n2 in itertools.combinations(list(graph.nodes()), 2):
        if n1[1] == n2[1]:
            graph.add_edge(n1, n2)
            expected_interactions.append((n1, n2))
    ordered_nodes = linear_graph_to_nodes_list(graph)
    positions = {
        node: (np.sin(a), np.cos(a))
        for (node, a) in zip(ordered_nodes,
                             np.linspace(0, 6.28, len(ordered_nodes) + 1))
    }


    graph = nx.Graph()

    for a in agents:
        graph.add_edge(*((a.name, site) for site in a.sites), rate=1.0,
                       is_slot=True, slot_name=a.name)
    for rule in rules:
        if rule.rate < rate_limit:
            continue
        graph.add_edge(*((rule.reactants[i].agent, rule.reactants[i].site)
                         for i in range(2)),
                       rate=rule.rate, is_interaction=True)
    
    if ax is None:
        fig, ax = plt.subplots(1, figsize=(len(slots), len(slots)))
    ax.axis("off")
    # ax.set_aspect('equal')
    for (n1, n2, data) in graph.edges(data=True):
        (x1, y1), (x2, y2) =  positions[n1], positions[n2]
        xmiddle, ymiddle = np.array((0.5 * (x1 + x2), 0.5 * (y1 + y2)))
        if data.get('is_interaction', False):
            if any(i in expected_interactions for i in [(n1, n2), (n2, n1)]):
                color, weight = 'grey', 'normal'
            else:
                color, weight = 'red', 'bold'
            ax.plot([x1, x2], [y1, y2], c=color, ls=':')
            factor = 1.2 if n1 == n2 else 1.0
            ax.text(factor * xmiddle, factor * ymiddle, "%d" % data['rate'],
                    ha='center', va='center',
                    fontdict=dict(color=color, weight=weight, size=10),
                    bbox=dict(boxstyle='round',
                          facecolor='#ffffffcc',
                          edgecolor='white'), zorder=1000)
        if data.get('is_slot', False):
            ax.plot([x1, x2], [y1, y2], c='b', ls='-', lw=6, alpha=0.3)
            ax.plot([xmiddle, 1.2 * xmiddle], [ymiddle, 1.3 * ymiddle],
                    c='b', ls='-', alpha=0.3, lw=2)
            ax.text(1.3 * xmiddle, 1.3 * ymiddle, data['slot_name'],
                    ha='center', va='center',
                    fontdict=dict(color='blue', weight='bold'),
                    bbox=dict(boxstyle='round',
                          facecolor='white',
                          edgecolor='white'), zorder=500)
            
            
    for node in graph.nodes():
        x, y = positions[node]
        ax.text(x, y, node[1], ha='center', va='center',
                bbox=dict(boxstyle='round',
                          facecolor='white',
                          edgecolor='white'))
    return ax