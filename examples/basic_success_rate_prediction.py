from kappagate import predict_assembly_accuracy, overhangs_list_to_slots


overhangs=  ['GGAG', 'GGCA', 'TCGC', 'CAGT', 'TCCA',
             'GAAT', 'AGTA', 'TCTT', 'CAAA', 'GCAC',
             'AACG', 'GTCT', 'CCAT']
slots = overhangs_list_to_slots(overhangs)
success_rate, _, _ = predict_assembly_accuracy(slots, initial_quantities=2000)
print ("Prediction for 12 High-Fi overhangs:", "%.1f%%" % (100 * success_rate))


overhangs = ['GGAG', 'GATA', 'GGCA', 'GGTC', 'TCGC',
             'GAGG', 'CAGT', 'GTAA', 'TCCA', 'CACA',
             'GAAT', 'ATAG', 'AGTA', 'ATCA', 'TCTT',
             'AGGT', 'CAAA', 'AAGC', 'GCAC', 'CAAC',
             'AACG', 'CGAA', 'GTCT', 'TCAG', 'CCAT']
slots = overhangs_list_to_slots(overhangs)
success_rate, _, _ = predict_assembly_accuracy(slots, initial_quantities=2000)
print ("Prediction for 24 High-Fi overhangs:", "%.1f%%" % (100 * success_rate))


overhangs=  ['GGAG', 'GGTC', 'AGCA', 'CAGT', 'GGTA',
             'GAAT', 'GGTT', 'TCTT', 'GGTG', 'GCAC',
             'AGCG', 'GTCT', 'CCAT']
slots = overhangs_list_to_slots(overhangs)
success_rate, _, _ = predict_assembly_accuracy(slots, initial_quantities=2000)
print ("Prediction for 12 Low-Fi overhangs:",
        "%.1f%%" % (100 * success_rate))
