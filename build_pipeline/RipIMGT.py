"""
A program to rip the sequences from genedb and parse them into fasta files

Ripped from here:

https://www.imgt.org/vquest/refseqh.html     
"""

from html.parser import HTMLParser
from html.entities import name2codepoint
import urllib.request, urllib.parse, urllib.error, os, sys


# Set globals
file_path  = os.path.split(__file__)[0]
html_outpath =  os.path.join( file_path, "IMGT_sequence_files", "htmlfiles" )
fasta_outpath = os.path.join( file_path, "IMGT_sequence_files", "fastafiles" )

# Define where to point the urls to.
# We have heavy, kappa, lambda, alpha, beta, gamma and delta chains.
# Both the v genes (imgt gapped amino acids) and the j genes (amino acids, are not gapped) 

# Urls as of 04-12-14
urls = { "HV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGHV&species=%s",
         "HJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGHJ&species=%s",
         "KV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGKV&species=%s",
         "KJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGKJ&species=%s",
         "LV": "https://www.imgt.org/genedb/GENElect?query=7.3+IGLV&species=%s",
         "LJ": "https://www.imgt.org/genedb/GENElect?query=7.6+IGLJ&species=%s",
         "AV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRAV&species=%s",
         "AJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRAJ&species=%s",
         "BV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRBV&species=%s",
         "BJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRBJ&species=%s",
         "GV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRGV&species=%s",
         "GJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRGJ&species=%s",
         "DV": "https://www.imgt.org/genedb/GENElect?query=7.3+TRDV&species=%s",
         "DJ": "https://www.imgt.org/genedb/GENElect?query=7.6+TRDJ&species=%s"
         
         #"HC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGHC&species=%s",
         #"KC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGKC&species=%s",
         #"LC": "https://www.imgt.org/genedb/GENElect?query=7.3+IGLC&species=%s",
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
           "Vicugna+pacos",
           "Bos+taurus"]


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
                print("1")
                try:
                    name = split[0]
                except IndexError:
                    return
                sequence = ("".join( split[1:])).replace(" ", "")
                self._data.append( (name, sequence) )

    def handle_entityref(self, name):
        self.currentnamedent = chr(name2codepoint[name])
    def handle_charref(self, name):
        if name.startswith('x'):
            self.currentnamedent = chr(int(name[1:], 16))
        else:
            self.currentnamedent = chr(int(name))

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
    if urllib.request.urlretrieve( urls[gene_type]%species,  filename ):
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
            print(">%s"%name, file=outfile)
            print(sequence, file=outfile)

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
            print("Bad parse", end=' ', file=sys.stderr)
            return 1
    else:
        print("Bad Url", end=' ', file=sys.stderr)
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
                print("Failed to retrieve %s %s"%(species, gene_type), file=sys.stderr)
            else:
                print("Parsed and saved %s %s"%(species, gene_type))
         
main()
