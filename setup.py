#!/usr/bin/env python3

#                               #
# Developed by James Dunbar     #
# Maintained by members of OPIG #
#                               #

from setuptools import setup

setup(name='anarci',
      version='1.3',
      description='Antibody Numbering and Receptor ClassIfication',
      author='James Dunbar',
      author_email='opig@stats.ox.ac.uk',
      url='http://opig.stats.ox.ac.uk/webapps/ANARCI',
      packages=['anarci'], 
      package_dir={'anarci': 'lib/python/anarci'},
      package_data={'anarci': ['dat/HMMs/ALL.hmm',
                               'dat/HMMs/ALL.hmm.h3f',
                               'dat/HMMs/ALL.hmm.h3i',
                               'dat/HMMs/ALL.hmm.h3m',
                               'dat/HMMs/ALL.hmm.h3p']},
      scripts=['bin/ANARCI']
     )

