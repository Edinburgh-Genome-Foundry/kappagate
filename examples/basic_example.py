from kappagate import predict_assembly_accuracy


overhangs=  ['GGAG', 'GGCA', 'TCGC', 'CAGT', 'TCCA', 'GAAT',
            'AGTA', 'TCTT', 'CAAA', 'GCAC', 'AACG', 'GTCT', 'CCAT']
predicted_rate, _, _ = predict_assembly_accuracy(overhangs=overhangs,
                                                 initial_quantities=2000)

print ("Prediction for 12 High-Fi overhangs:",
       "%.1f%%" % (100 * predicted_rate))



overhangs = ['GGAG', 'GATA', 'GGCA', 'GGTC', 'TCGC', 'GAGG',
             'CAGT', 'GTAA', 'TCCA', 'CACA', 'GAAT', 'ATAG',
             'AGTA', 'ATCA', 'TCTT', 'AGGT', 'CAAA', 'AAGC',
             'GCAC', 'CAAC', 'AACG', 'CGAA', 'GTCT', 'TCAG', 'CCAT']
predicted_rate, _, _ = predict_assembly_accuracy(overhangs=overhangs,
                                                 initial_quantities=2000)

print ("Prediction for 24 High-Fi overhangs:",
       "%.1f%%" % (100 * predicted_rate))



overhangs=  ['GGAG', 'GGTC', 'AGCA', 'CAGT', 'GGTA', 'GAAT', 'GGTT',
             'TCTT', 'GGTG', 'GCAC', 'AGCG', 'GTCT', 'CCAT']
predicted_rate, _, _ = predict_assembly_accuracy(overhangs=overhangs,
                                                 initial_quantities=2000)

print ("Prediction for 12 Low-Fi overhangs:",
        "%.1f%%" % (100 * predicted_rate))
