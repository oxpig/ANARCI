"""
A program to rip the sequences from genedb and parse them into fasta files

Ripped from here:

http://www.imgt.org/vquest/refseqh.html     
"""

from HTMLParser import HTMLParser
from htmlentitydefs import name2codepoint
import urllib, os, sys


# Set globals
file_path  = os.path.split(__file__)[0]
html_outpath =  os.path.join( file_path, "IMGT_sequence_files", "htmlfiles" )
fasta_outpath = os.path.join( file_path, "IMGT_sequence_files", "fastafiles" )

# Define where to point the urls to.
# We have heavy, kappa, lambda, alpha, beta, gamma and delta chains.
# Both the v genes (imgt gapped amino acids) and the j genes (amino acids, are not gapped) 

# Urls as of 04-12-14
urls = { "HV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGHV&species=%s",
         "HJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+IGHJ&species=%s",
         "KV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGKV&species=%s",
         "KJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+IGKJ&species=%s",
         "LV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGLV&species=%s",
         "LJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+IGLJ&species=%s",
         "AV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRAV&species=%s",
         "AJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRAJ&species=%s",
         "BV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRBV&species=%s",
         "BJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRBJ&species=%s",
         "GV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRGV&species=%s",
         "GJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRGJ&species=%s",
         "DV": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+TRDV&species=%s",
         "DJ": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.6+TRDJ&species=%s"
         
         #"HC": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGHC&species=%s",
         #"KC": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGKC&species=%s",
         #"LC": "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.3+IGLC&species=%s",
       }



# Species as of 04-12-14
# Species as of 02-06-16 - alpaca added
# These are retrieved for all the antibodies
all_species = ["Homo+sapiens",
           "Mus",
           "Rattus+norvegicus",
           "Oryctolagus+cuniculus",
           "Macaca+mulatta",
           "Sus+scrofa",
           "Vicugna+pacos",]


# These are retrieved for the tcr chains. There are a few more for gamma and delta chains but
# they are rare (structurally anyway) that it does not seem worth it.   
all_tr_species = ["Homo+sapiens",
           "Mus",
]


# These do not have light chain sequences so we ignore (they're fish)
#           "Oncorhynchus+mykiss",
#           "Danio+rerio" ]

# Html parser class.
class GENEDBParser(HTMLParser):
    currenttag = None
    currentnamedent = None
    _data = []
    def handle_starttag(self, tag, attrs):
        self.currenttag=tag
    def handle_endtag(self, tag):
        self.currenttag=None
    def handle_data(self, data):
        split = data.split("\n")
        start = sum([ 1 if l[0]==">" else 0 for l in split if len(l)])
        if self.currenttag=="pre" and (self.currentnamedent ==">" or start):
            # Two different ways of parsing the html based on how IMGT have formatted the pages.
            # For some reason they format gene db differently sometimes (legacy?) 
            if start > 1: # If you encounter more than one line in the data with a fasta ">" symbol, all sequences will be in the same packet
                name, sequence = None, ""
                for l in split:
                    if not l: continue
                    if l[0]==">":
                        if sequence:
                            self._data.append( (name, sequence) )
                            name, sequence = None, ""
                        name = l
                    else:
                        sequence += l.replace(" ", "")
                if name and sequence:
                    self._data.append( (name, sequence) )
            else: # Otherwise it will be done entry by entry
                print "1"
                try:
                    name = split[0]
                except IndexError:
                    return
                sequence = ("".join( split[1:])).replace(" ", "")
                self._data.append( (name, sequence) )

    def handle_entityref(self, name):
        self.currentnamedent = unichr(name2codepoint[name])
    def handle_charref(self, name):
        if name.startswith('x'):
            self.currentnamedent = unichr(int(name[1:], 16))
        else:
            self.currentnamedent = unichr(int(name))

    def rip_sequences(self,htmlstring):
        """
        Method for this subclass that automates the return of data
        """
        self.reset()
        self._data = []
        self.currenttag = None
        self.currentnamedent = None
        self.feed(htmlstring)
        self._data
        return self._data
        
parser = GENEDBParser()


def get_html(species, gene_type, force = True):
    """
    Get the html from IMGT
    """
    filename = os.path.join(html_outpath,"%s_%s.html"%(species.replace("+", "_"), gene_type) )
    # If html file exists already
    if os.path.isfile(filename):
        return filename
    if urllib.urlretrieve( urls[gene_type]%species,  filename ):
        return filename
    else:
        return False

def write_fasta( sequences, species, gene_type ):
    """
    Write a fasta file containing all sequences
    """
    filename = os.path.join(fasta_outpath,"%s_%s.fasta"%(species.replace("+", "_"), gene_type) )
    with open(filename, "w") as outfile:
        for name, sequence in sequences:
            print >> outfile, ">%s"%name
            print >> outfile, sequence

def ripfasta(species, gene_type, force = True):
    """ 
    Rip the fasta sequences for a species and gene type from IMGT
    """
    htmlfile = get_html(species, gene_type, force)
    if htmlfile:
        with open(htmlfile) as infile:
            sequences = parser.rip_sequences(infile.read())
        if sequences:
            write_fasta(sequences, species, gene_type )
        else:
            print >> sys.stderr, "Bad parse",
            return 1
    else:
        print >> sys.stderr, "Bad Url",
        return 1

def main():
    """
    For all V and J gene types (H,K,L,A,B,G,D) parse IMGT database and extract fasta files
    """
    for gene_type in urls:
        for species in all_species:
            if gene_type[0] in "ABGD" and species not in all_tr_species: continue # we don't want TCRs for all organisms
            if gene_type[0] in "KL" and species == "Vicugna+pacos": continue # alpacas don't have light chains
            if ripfasta(species, gene_type, force = False):
                print >> sys.stderr, "Failed to retrieve %s %s"%(species, gene_type)
            else:
                print "Parsed and saved %s %s"%(species, gene_type)
         
main()
