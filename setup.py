import shutil, os, subprocess
import site, sys
from importlib import util
from distutils.core import setup
from setuptools.command.install import install

class CustomInstallCommand(install):
   def run(self):
       install.run(self)
       # Post-installation routine
       ANARCI_LOC = os.path.join(site.getsitepackages()[0], 'anarci') # site-packages/ folder
       ANARCI_BIN = sys.executable.split('python')[0] # bin/ folder

       shutil.copy('bin/ANARCI', ANARCI_BIN) # copy ANARCI executable
       print("INFO: ANARCI lives in: ", ANARCI_LOC) 

       # Build HMMs from IMGT germlines
       os.chdir("build_pipeline")
       print('INFO: Downloading germlines from IMGT and building HMMs...')
       print("INFO: running 'RUN_pipeline.sh', this will take a couple a minutes.")
       proc = subprocess.Popen(["bash", "RUN_pipeline.sh"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
       o, e = proc.communicate()

       print(o.decode())
       print(e.decode())
       
       # Copy HMMs where ANARCI can find them
       shutil.copy( "curated_alignments/germlines.py", ANARCI_LOC )
       os.mkdir(os.path.join(ANARCI_LOC, "dat"))
       shutil.copytree( "HMMs", os.path.join(ANARCI_LOC, "dat/HMMs/") )
      
      # Remove data from HMMs generation
       try:
           shutil.rmtree("curated_alignments/")
           shutil.rmtree("muscle_alignments/")
           shutil.rmtree("HMMs/")
           shutil.rmtree("IMGT_sequence_files/")
       except OSError:
           pass

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
     scripts=['bin/ANARCI'],
     cmdclass={"install": CustomInstallCommand, }, # Run post-installation routine
    )
