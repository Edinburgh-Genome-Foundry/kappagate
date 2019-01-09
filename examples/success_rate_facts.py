from kappagate import (overhangs_list_to_slots, predict_assembly_accuracy,
                       plot_colony_picking_graph, success_rate_facts)

overhangs = ['TAGG', 'GACT', 'GGAC', 'CAGC',
             'GGTC', 'GCGT', 'TGCT', 'GGTA',
             'CGTC', 'CTAC', 'GCAA', 'CCCT']
slots = overhangs_list_to_slots(overhangs)
predicted_rate, _, _ = predict_assembly_accuracy(slots)
ax = plot_colony_picking_graph(success_rate=predicted_rate)
ax.figure.savefig("success_rate_facts.png", bbox_inches='tight')

print("SUMMARY:")
print (success_rate_facts(predicted_rate, plain_text=True))
print ("See 'success_rate_facts.png' for an illustration.")
