from kappagate import overhangs_list_to_slots, plot_circular_interactions
overhangs = ['TAGG', 'GACT', 'GGAC', 'CAGC',
             'GGTC', 'GCGT', 'TGCT', 'GGTA',
             'CGTC', 'CTAC', 'GCAA', 'CCCT']
slots = overhangs_list_to_slots(overhangs)
ax = plot_circular_interactions(
    slots, annealing_data=('25C', '01h'), rate_limit=200)
ax.figure.savefig("plotting_interactions.png", bbox_inches='tight')