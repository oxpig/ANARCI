#! /usr/bin/env python

import shutil, os, subprocess, importlib
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
      packages=['anarci'], 
      package_dir={'anarci': 'lib/python/anarci'},
      data_files = [ ('bin', ['bin/muscle', 'bin/muscle_macOS', 'bin/ANARCI']) ],
      include_package_data = True,
      scripts=['bin/ANARCI']
     )


ANARCI_LOC = importlib.util.find_spec("anarci").submodule_search_locations[0]
os.chdir("build_pipeline")

print('Downloading germlines from IMGT and building HMMs...')
proc = subprocess.Popen(["bash", "RUN_pipeline.sh"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
o, e = proc.communicate()

print(o.decode())
print(e.decode())

shutil.copy( "curated_alignments/germlines.py", ANARCI_LOC )
shutil.copytree( "HMMs", os.path.join(ANARCI_LOC, "dat/HMMs/") )

try:
    shutil.rmtree("curated_alignments/")
    shutil.rmtree("muscle_alignments/")
    shutil.rmtree("HMMs/")
    shutil.rmtree("IMGT_sequence_files/")
    os.mkdir(os.path.join(ANARCI_LOC, "dat"))
except OSError:
    pass

