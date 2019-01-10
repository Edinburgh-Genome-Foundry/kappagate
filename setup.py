import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('kappagate/version.py').read()) # loads __version__

setup(
    name='kappagate',
    version=__version__,
    author='Zulko',
    url='https://github.com/Edinburgh-Genome-Foundry/kappagate',
    description='Predict valid-clone-rates in Golden Gate DNA assemblies',
    long_description=open('pypi-readme.rst').read(),
    license='MIT',
    keywords="DNA assembly synthetic biology golden gate",
    packages=find_packages(exclude='docs'),
    install_requires=['topkappy', 'networkx', 'tatapov', 'matplotlib',
                      'dnacauldron', 'proglog', 'flametree', 'biopython',
                      'snapgene_reader'])
