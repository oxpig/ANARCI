#!/usr/bin/env python

#                               #
# Developed by James Dunbar     #
# Maintained by members of OPIG #
#                               #

import shutil, os, subprocess, imp
# Clean this out if it exists
if os.path.isdir("build"):
    shutil.rmtree("build/")

from distutils.core import setup

setup(name='anarci',
      version='1.3',
      description='Antibody Numbering and Receptor ClassIfication',
      author='James Dunbar',
      author_email='opig@stats.ox.ac.uk',
      url='http://opig.stats.ox.ac.uk/webapps/ANARCI',
      packages=['anarci', 'anarci.Bio', 'anarci.Bio.Align',
                #'anarci.builder',
      'anarci.Bio.Alphabet', 'anarci.Bio.Data', 'anarci.Bio._py3k',
      'anarci.Bio.SearchIO','anarci.Bio.SearchIO._model', 'anarci.Bio.SeqIO', 'anarci.Bio.SearchIO.HmmerIO', 'anarci.Bio.AlignIO'  ],
      package_dir={'anarci': 'lib/python/anarci',
                 #  'anarci.builder': 'lib/python/anarci/builder',
                   'anarci.Bio.AlignIO': 'lib/python/anarci/Bio/AlignIO', 
                   'anarci.Bio.SearchIO': 'lib/python/anarci/Bio/SearchIO', 
                   'anarci.Bio.SearchIO.HmmerIO': 'lib/python/anarci/Bio/SearchIO/HmmerIO', 
                   'anarci.Bio.Alphabet': 'lib/python/anarci/Bio/Alphabet', 
                   'anarci.Bio._py3k': 'lib/python/anarci/Bio/_py3k', 
                   'anarci.Bio.Data': 'lib/python/anarci/Bio/Data', 
                   'anarci.Bio': 'lib/python/anarci/Bio', 
                   'anarci.Bio.SearchIO._model': 'lib/python/anarci/Bio/SearchIO/_model', 
                   'anarci.Bio.SeqIO': 'lib/python/anarci/Bio/SeqIO',
                   'anarci.Bio.Align': 'lib/python/anarci/Bio/Align'},
      #package_data={'anarci': ['dat/HMMs/ALL.hmm',
      #                         'dat/HMMs/ALL.hmm.h3f',
      #                         'dat/HMMs/ALL.hmm.h3i',
      #                         'dat/HMMs/ALL.hmm.h3m',
      #                         'dat/HMMs/ALL.hmm.h3p']},
      scripts=['bin/ANARCI', "bin/muscle", "bin/muscle_macOS"],
      license="GPLv3"
     )

####
import sys
if sys.argv[1] != "install":
    sys.exit(0)

try:
    ANARCI_LOC = imp.find_module("anarci")[1]
except:
    sys.stderr.write("Something isn't right. Aborting.")
    sys.exit(1)

os.chdir("build_pipeline")

try:
    shutil.rmtree("curated_alignments/")
    shutil.rmtree("muscle_alignments/")
    shutil.rmtree("HMMs/")
    shutil.rmtree("IMGT_sequence_files/")
    os.mkdir(os.path.join(ANARCI_LOC, "dat"))
except OSError:
    pass

print 'Downloading germlines from IMGT and building HMMs...'
proc = subprocess.Popen(["bash", "RUN_pipeline.sh"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
o, e = proc.communicate()

print o
print e

shutil.copy( "curated_alignments/germlines.py", ANARCI_LOC )
shutil.copytree( "HMMs", os.path.join(ANARCI_LOC, "dat/HMMs/") )
