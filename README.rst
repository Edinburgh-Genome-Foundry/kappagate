.. raw:: html

    <p align="center">
    <img alt="kappagate Logo" title="kappagate Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/kappagate/master/docs/title.png" width="600">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/kappagate.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/kappagate
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/kappagate/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/kappagate?branch=master


Kappagate is a Python library to predict the percentage of good clones (carrying
a correct version of the desired assembly) when assembling DNA with a method
relying on 4bp overhangs (e.g. Golden Gate assembly, OGAB, etc.).

Using Kappagate, you can get an estimation of how difficult the
assembly will be, and how many clones should be tested to find a correct one. 

Kappagate uses the exhaustive relative overhangs affinities tables provided
in Potapov et. al. 2018 (ACS Syn. Bio.). In this publication the authors show
that the proportion of valid clones rates can be predicted using focused
in-vitro experiments focused on the overhangs present in the assembly.

Kappagate attempts to predict clone validity rates without any overhangs-subset-specific
experiment, using computer simulations instead. It simulates the temporal evolution
of the DNA fragments ligation reaction using the Kappa biological modeling system.
At the end of the cloning simulation, Kappagate returns the ratio between "good"
constructs (with all expected parts in the right order) and bad circular assembly-forming
constructs (which may produce bad clones after transformation and plating).

This is an experimental piece of software, useful to us, but coming with no warranty.

Example
--------


.. code:: python

    from kappagate import predict_assembly_accuracy

    # FIRST TEST ON 12 WELL-DESIGNED OVERHANGS

    overhangs=  ['GGAG', 'GGCA', 'TCGC', 'CAGT', 'TCCA', 'GAAT',
                 'AGTA', 'TCTT', 'CAAA', 'GCAC', 'AACG', 'GTCT', 'CCAT']
    predicted_rate, _, _ = predict_assembly_accuracy(overhangs=overhangs)
    
    print (predicted_rate)
    # >>> 0.987

This means that 98.7% of clones will carry a valid assembly. It is really
not far from the experimental observation in,, Potapov et al. which was 99.2% +- 0.6.
Let's have a look at a few more sets:

.. code:: python

    overhangs = ['GGAG', 'GATA', 'GGCA', 'GGTC', 'TCGC', 'GAGG',
                 'CAGT', 'GTAA', 'TCCA', 'CACA', 'GAAT', 'ATAG',
                 'AGTA', 'ATCA', 'TCTT', 'AGGT', 'CAAA', 'AAGC',
                 'GCAC', 'CAAC', 'AACG', 'CGAA', 'GTCT', 'TCAG', 'CCAT']
    predicted_rate, _, _ = predict_assembly_accuracy(overhangs=overhangs)
    print (predicted_rate)
    # >>> 0.84
    # In Potapov 2018: 84% +/- 5%

.. code:: python

    overhangs=  ['GGAG', 'GGTC', 'AGCA', 'CAGT', 'GGTA', 'GAAT', 'GGTT',
                'TCTT', 'GGTG', 'GCAC', 'AGCG', 'GTCT', 'CCAT']
    predicted_rate, _, _ = predict_assembly_accuracy(overhangs=overhangs)
    print (predicted_rate)
    # >>> 0.33
    # In Potapov 2018: 45% +/- 5%

Installation
-------------

You can install kappagate through PIP

.. code::

    sudo pip install kappagate

Alternatively, you can unzip the sources in a folder and type

.. code::

    sudo python setup.py install

License = MIT
--------------

Kappagate is an open-source software originally written at the `Edinburgh Genome Foundry <http://genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_ and `released on Github <https://github.com/Edinburgh-Genome-Foundry/kappagate>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

More biology software
---------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

Kappagate is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
