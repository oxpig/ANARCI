#!/usr/bin/env python3

#                               #
# Developed by James Dunbar     #
# Maintained by members of OPIG #
#                               #

import shutil, os, subprocess, importlib

# Clean this out if it exists
if os.path.isdir("build"):
    shutil.rmtree("build/")

from distutils.core import setup
from setuptools.command.install import install


def build():
    try:
        ANARCI_LOC = os.path.dirname(importlib.util.find_spec("anarci").origin)
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

    print('Downloading germlines from IMGT and building HMMs...')
    proc = subprocess.Popen(["bash", "RUN_pipeline.sh"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    o, e = proc.communicate()

    print(o.decode())
    print(e.decode())

    shutil.copy( "curated_alignments/germlines.py", ANARCI_LOC )
    shutil.copytree( "HMMs", os.path.join(ANARCI_LOC, "dat/HMMs/") )


class Install(install):

    def run(self):
        super().run()
        build()


setup(name='anarci',
      version='1.3',
      description='Antibody Numbering and Receptor ClassIfication',
      author='James Dunbar',
      author_email='opig@stats.ox.ac.uk',
      url='http://opig.stats.ox.ac.uk/webapps/ANARCI',
      packages=['anarci'], 
      package_dir={'anarci': 'lib/python/anarci'},
      #package_data={'anarci': ['dat/HMMs/ALL.hmm',
      #                         'dat/HMMs/ALL.hmm.h3f',
      #                         'dat/HMMs/ALL.hmm.h3i',
      #                         'dat/HMMs/ALL.hmm.h3m',
      #                         'dat/HMMs/ALL.hmm.h3p']},
      scripts=['bin/ANARCI'],
      install_requires=['biopython>=1.78'],
      data_files = [ ('bin', ['bin/muscle', 'bin/muscle_macOS']) ],
      cmdclass={'install': Install}
     )


