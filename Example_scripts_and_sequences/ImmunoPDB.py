#! /usr/bin/env python
#    ANARCI - Antibody Numbering and Antigen Receptor ClassIfication
#    Copyright (C) 2016 Oxford Protein Informatics Group (OPIG)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the BSD 3-Clause License.
#
#    You should have received a copy of the BSD 3-Clause Licence
#    along with this program.  If not, see <https://opensource.org/license/bsd-3-clause/>.

description='''

ANARCI - ImmunoPDB                                     \\\    //
Antibody Numbering and Antigen Receptor ClassIfication  \\\  //
                                                          ||
(c) Oxford Protein Informatics Group (OPIG). 2015-16      ||

Example script to number Antibody and TCR PDB structures.

ANARCI must be installed and in the path
opig.stats.ox.ac.uk/webapps/anarci

Requirements: Biopython (version >= 1.66)
              Muscle


This script extends the BioPython PDBParser and Structure classes so that a numbering scheme can be applied to the variable domain
of an antigen receptor chain. 

o Where available the Seqres record is used as the full sequence. Missing residues are recognised by comparing this to the residues
  with coordinates.

o *Only* variable domains are numbered consistently in the chosen scheme.

o Residues before the domain are numbered '0' with reverse alphabetical insertions if there are less than 28 (all 0 otherwise - this
  will break some PDB parsers...)

o Residues after the variable domain are numbered sequentially from 1001.

o By default when more than one variable domain is found on a single chain (e.g. single chain Fv, diabody...) the numbering will be
  with respect to the first domain identified.

o CDR recognition is performed and regions are annotated in the xtra dictionary attributes of residue objects.

o Pairing is performed using the distance between the interface cysteine positions (imgt 104).

Basic useage

Renumber antibody chains with imgt numbering scheme
python ImmunoPDB.py -i infile.pdb -o outfile.pdb -s imgt

Renumber tcr chains with imgt numbering scheme
python ImmunoPDB.py -i infile.pdb -o outfile.pdb -s imgt --receptor tr


'''

epilogue='''
Author: James Dunbar (dunbar@stats.ox.ac.uk)
        Charlotte Deane (deane@stats.ox.ac.uk)

Copyright (C) 2016 Oxford Protein Informatics Group (OPIG)
Freely distributed under the BSD 3-Clause Licence.
'''

# Antigen receptor numbering and classification.
from anarci import anarci, scheme_names, scheme_short_to_long

# Biopython - you must use the dev version for this to work. Got to: https://github.com/biopython/biopython
# Must be version 1.84 or greater.
from Bio.PDB import *
from Bio.File import as_handle
# from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SCOP.Raf import protein_letters_3to1_extended as protein_letters_3to1
from Bio.SeqUtils import seq1

# Python
import subprocess
import os
import collections


# Globals
tebahpla = 'ZYXWVUTSRQPONMLKJIHGFEDCBA '
cysh = {'kabat':92,'chothia':92,'martin':92,'imgt':104,'aho':106,'wolfguy':330 } # The position of the interface cysteine in each
cysl = {'kabat':88,'chothia':88,'martin':88,'imgt':104,'aho':106,'wolfguy':734 } # scheme.
inter_cys = { 'H':cysh, 'B':cysh, 'D':cysh, 'L':cysl, 'K':cysl, 'A':cysl, 'G':cysl }


# Class to apply numbering to a structure object. This can be used in conjunction with the Ig/Tcr parsers or applied to an existing
# object. This class annotates in the .xtra dictionary of chains and residues and *does not* modify the residue identifiers.
class PDBNumber(object):
    '''
    This will annotate residues in the structure with numbering schemes
    The original structure object remains unchanged.
    Annotations are listed in the residue .xtra dictionary. 
    Pairings are listed in the Chains .xtra dictionary.
    '''
    def __init__(self,scheme,allowed_domains, region_definition=None, warnings=False):
        self.scheme=scheme_short_to_long[scheme.lower()]
        self.allowed_domains = allowed_domains
        self.region_definition = region_definition
        self.warnings = warnings

    def numberChain(self, chain, sequence=''):
        '''
        Number a Bio.PDB chain with a Ig or Tcr numbering scheme.
        '''
        check_continuity = False
        structureSequence, index2ResId = extract_sequence( chain )
        if not sequence: # No seqres provided
            sequence = structureSequence
            seqI2strucI = get_alignment_dict( sequence, sequence )
            if self.warnings: check_continuity = True
        elif which( 'muscle' ) or structureSequence in sequence: # Muscle found or it is an easy alignment
            chain.xtra['structure_ali'], chain.xtra['seqres_ali'] = pairwise_muscle( structureSequence , sequence ) 
            fix_deletions( chain ) # Fix the position of a missing residue if a degenerate alignment
            seqI2strucI = get_alignment_dict(chain.xtra['seqres_ali'], chain.xtra['structure_ali'])
        else: # Missing residues are suspected. 
            raise Warning('Muscle was not found and there are suspected missing residues. Numbering could be incorrect')
            sequence = structureSequence
            seqI2strucI = get_alignment_dict( sequence, sequence )

        # Reset if needed
        if 'anarci_details' in chain.xtra:
            chain.xtra['domains'] = []
            chain.xtra['anarci_details'] = []

        di = 0 # Domain index (when there are multiple variable domains on the same chain)
        for numbering2Index, domainType, details in self.numberSequence(chain.id, sequence): # Iterate over the recognised domains.
            if domainType in self.allowed_domains:
                chain.xtra['domains'] = chain.xtra.get('domains', [] ) + [domainType] # Append to the domains list.
                chain.xtra['anarci_details'] = chain.xtra.get('anarci_details', [] ) + [details] # Append to the germlines list.
                for numberedId, index in numbering2Index:
                    if index not in seqI2strucI: continue # Missing residue
                    res = chain[ index2ResId[ seqI2strucI[index] ] ]
                    res.xtra[self.scheme] = numberedId # Store the annotated numbering               
                    if 'pdb' not in res.xtra: res.xtra['pdb']  = res.id[1:] # Store the original numbering
                    res.xtra['domain_type'] = domainType # Store the domain type
                    res.xtra['domain_index'] = di # Store which domain this residue is in.

                # Second iteration to annotate non variable domains.
                if 1 < details['query_start'] < 28: # If possible, annotate pre domain with 0C 0B 0A as pre insertions. 
                    insertions = tebahpla[ 27 - details['query_start']: ]
                else:
                    insertions = ' '*(details['query_start']+1) # So that insertion code will always be at least ' ' # All 
                
                pre, ri, ii  = True, 0, 0
                for residue in chain:
                    if 'pdb' not in residue.xtra: 
                        residue.xtra['pdb'] = residue.id[1:]
                    if residue.xtra.get('domain_index','') == di:
                        if pre:
                            pre, ri, ii = False, 1001, 0
                        continue
                    residue.xtra[self.scheme+'%d'%di] = (ri, insertions[ii])
                    if pre: # Update the insertion code 
                        ii += 1
                    else: # Update the residue index 
                        ri += 1                    
                self.annotateChainRegions( chain, self.scheme, self.region_definition, domain_index=di )
                di +=1

        if di and check_continuity:
            continuous, at = analyse_continuity( chain )
            if not continuous:
                print('Warning: Numbering may be incorrect as a missing residue was detected in chain %s (around residue %d%s). Provide the seqres record to overcome this problem.'%(chain.id, chain.child_list[at].id[1], chain.child_list[at].id[2] ))
            
        chain.xtra['scheme'] = self.scheme

    def annotateChainRegions(self, chain, scheme, definition=None, domain_index=0):
        '''
        Annotate chain with regions. The loaded numbering scheme is used
        '''
        if definition is None or scheme=='pdb': return
    
        domain_type = chain.xtra['domains'][domain_index].replace('K','L')
        numbering   = [ (r.xtra[scheme], None) if r.xtra.get('domain_index',None) == domain_index else ((0,' '), None) for r in chain ]

        try:
            regions = annotate_regions(numbering, domain_type, numbering_scheme=scheme,definition=definition)
        except AssertionError: # Unimplemented
            return
        for i in range(len(numbering)):
            if not chain.child_list[i].xtra.get('region'): # Update the region.
                chain.child_list[i].xtra['region'] = regions[i][2]

    def numberSequence(self,sequenceName, sequence):
        '''
        Number a sequence and yield the annotations for each domain identified (n to c)
        '''
        if len(sequence) > 70: 
            results = anarci( [(sequenceName, sequence)], scheme=self.scheme, assign_germline=True, allow=self.allowed_domains )
        else:
            return 
        numbered, details = results[0][0], results[1][0]
        if numbered is None: numbered = []
        for i in range( len( numbered ) ): # Iterate over the identified domains (e.g. for an scfv)
            numbering = [ (n, a) for n, a in numbered[i][0] if a != '-' ] # Remove gaps if made (imgt scheme)
            yield [ (numbering[ri][0], ri+numbered[i][1]) for ri in range( len( numbering ) ) ], details[i]['chain_type'], details[i]

# Extension of the Biopython PDBParser. 
class AntigenReceptorPDBParser(PDBParser):
    '''
    A generic class to parse antigen receptor structures. 
    '''
    allowed_schemes = scheme_names
    allowed_domains = set(['H','K','L'] + ['A','B','G','D'])
    allowed_pairs = [ ('H','K'),('H','L'),('B','A'),('D','G') ]
    annotates = 'antigen receptor'        
    
    def __init__( self, PERMISSIVE=True, get_header=False,
                 structure_builder=None, QUIET=False, scheme='imgt', region_definition='imgt', warnings=False):
        '''
        '''
        super( AntigenReceptorPDBParser, self ).__init__(PERMISSIVE=PERMISSIVE, get_header=get_header,
                 structure_builder=structure_builder, QUIET=QUIET)
        assert scheme in self.allowed_schemes, 'Scheme unknown'
        self.scheme=scheme_short_to_long[scheme.lower()]
        self.region_definition = region_definition
        self.warnings = warnings

    def get_structure(self, id, file):
        """Return the structure. Variable domains are annotated with the chosen numbering scheme.

        o Receptor chains are numbered with the desired numbering scheme and  paired according to the separation of their interface
        cysteine (or whatever is at position 104 in imgt). 

        o Only 'opposite' type domains (VH-VK, VH-VL, VB-VA, VD-VG) are paired. VL-VL dimers will not be annotated

        o Single chain fvs are also considered.

        o The xtra['pairing'] list of a chain tells you how the variable domains it contains are paired. 
            e.g. For a normal Fv with chains H and L. 
                structure[0]['H'].xtra['pairing']  -->    [ (<Chain id=L>, 0) ]
                structure[0]['L'].xtra['pairing']  -->    [ (<Chain id=H>, 0) ]
                                                         chain obj,  domain index of paired              

            e.g. For a scfv with H and L domains on the same chain, A
                structure[0]['A'].xtra['pairing']   -->   [ (<Chain id=A>, 1) , (<Chain id=A>, 0) ]

        o Numbering will be returned with respect to the first identified variable domain. Use the domain_index argument in the 
          structure object's 'switch_numbering_scheme' method to change to number with respect to the second etc for an scfv.
        
        o Residues before the variable domain are all given the residue number 0. If there are less than 27 residues before the start
          of the variable domain they will be annotated with insertion codes in reverse alphabetical order
            e.g.
             0D 0C 0B 0A 0 | 1 2 3 4
        
        o Residues after the variable domain are numbered sequentially from 1001. No guarantee is given of structural equivalence in 
          the constant domain (yet...) 


        Arguments:
         - id - string, the id that will be used for the structure
         - file - name of the PDB file OR an open filehandle
        """

        # Ensure the file is open file handle 
        with as_handle( file ) as fin:

            # Read the structure using biopyton
            structure = super( AntigenReceptorPDBParser, self ).get_structure(id, fin)

            # Rewind to the start of the file
            fin.seek(0)

            # Read the seqres entries
            seqres    = self.read_seqres( fin )

        stringseqres = dict( (e, str(seqres[e].seq) ) for e in seqres )
        self.number_receptor_chains(structure, sequences=stringseqres) # Post process
        for chain in seqres: structure[0][chain].xtra['seqresobj'] = seqres[chain] # Associtate with the chain object - includes the original lines
        self.find_pairs(structure) # Pair chains
        structure.switch_numbering_scheme(self.scheme) # Return the structure with the requested numbering scheme. 
        return structure

    def read_seqres(self, file):
        '''
        Read the seqres entries of a structure if available
        '''
        try:
            with as_handle(file) as f:
                seqres   =  dict([ (a.id, a) for a in PdbSeqresIterator(f) ])
        except:
            seqres = {}
        return seqres

    def number_receptor_chains(self, structure, sequences={}):
        '''
        Apply the numbering to each PDB chain in the file
        '''
        pdbN = PDBNumber(self.scheme, self.allowed_domains, region_definition=self.region_definition, warnings=self.warnings)
        for model in structure:
            for chain in model:
                pdbN.numberChain(chain, sequence=sequences.get(chain.id, '') ) 
                chain.xtra['seqres'] = sequences.get(chain.id, '')


    def find_pairs(self,structure): 
        '''
        Annotate the pairs of domains in the structure using the distance between the interface cysteine residue CA atoms.
        '''
        for model in structure:
            # Bin the domains of the structure by type
            model.xtra['domains'] = dict( (d, []) for d in self.allowed_domains )
            model.xtra['firstpair']=None # Used for reference if outputting only the first receptor identified.
            for chain in model:
                for di in range( len( chain.xtra.get('domains',[]) ) ):
                    model.xtra['domains'][ chain.xtra['domains'][di] ].append( (chain, di) )
                    if not model.xtra['firstpair']: model.xtra['firstpair'] = [(chain, di)] # For single domain
                chain.xtra['pairing'] = [None]*len( chain.xtra.get('domains',[]) )
            for pA, pB in self.allowed_pairs:
               
               for chainA, diA in model.xtra['domains'][pA]:
                    aA = chainA.interface_cys_CA( self.scheme, diA )
                    if aA is None: continue
                    for chainB, diB in model.xtra['domains'][pB]:
                        aB = chainB.interface_cys_CA( self.scheme, diB )
                        if aB is None: continue
                        if aA - aB < 25: 
                            chainA.xtra['pairing'][diA] = (chainB, diB) 
                            chainB.xtra['pairing'][diB] = (chainA, diA)
                            if len( model.xtra['firstpair']) < 2: model.xtra['firstpair'] = sorted( [(chainA, diA) ,(chainB, diB)] )

class TcrPDBParser(AntigenReceptorPDBParser):
    allowed_schemes = ['imgt','i','aho','a']
    allowed_domains  = set(['A','B','G','D'])
    allowed_pairs = [ ('B','A'),('D','G') ] 
    annotates = 'tcr'

class AntibodyPDBParser(AntigenReceptorPDBParser):
    allowed_domains  = set(['H','K','L'])
    allowed_pairs = [ ('H','K'),('H','L')]
    annotates = 'antibody'

# Only accept variable regions of a chain
class SelectFv(Select):
    def accept_residue(self, residue):
        if 'anarci_details' in residue.parent.xtra:
            if residue.id[1] > 0 and residue.id[1] < 1000:
                return 1

# Allows non-receptor chains but restricts scfvs to the variable part.   
class SelectFvScFv(Select):
    def accept_residue(self, residue):
        if len(residue.parent.xtra.get('domains',[])) > 1:
            if residue.id[1] > 0 and residue.id[1] < 1000:
                return 1
            return 0 
        return 1


# Extend the Biopython classes with methods to handle numbering.
def Entity_switch_numbering_scheme(self, scheme, domain_index=0): # Structure, Model
    ''' 
    '''
    for child in self:
        child.switch_numbering_scheme(scheme, domain_index=domain_index)

def Chain_switch_numbering_scheme(self, scheme, domain_index=0):
    '''
    Change the numbering scheme used in the identifiers of the residues in the chain. 
    By default the numbering will be with respect to the first (0th) variable domain identified.
    If a single chain fv is identified (i.e. two variable domains on the same chain) one can switch the numbering to be with
    respect to the second domain using the domain index argument (0 - first domain, 1 - second domain etc). 
    
    Residues before the numbered domain are numbered 0. Insertion codes are used in reverse alphabetical order iff there are 27 or
    less n-terminal 'insertions' (all 0 otherwise). Residues after the domain are numbered sequentially from 1001.
   
    '''
    if scheme not in [ self.xtra['scheme'], 'pdb' ]:
        PDBNumber(scheme, self.xtra.get('domains', ['H','K','L','A','B','G','D'])).numberChain( self, sequence=self.xtra.get('seqres','' ) )

    if 'anarci_details' in self.xtra: # Only bother to recurse into numbered chains 
        if self.xtra.get( 'loaded','' ) == scheme+str(domain_index): return # Only switch if different.
        for child in self:
            child.switch_numbering_scheme(scheme, domain_index=domain_index)
        self.child_dict = {} # Reset the child dict on a separate loop to prevent conflicts
        for child in self:  
            self.child_dict[child.id]=child
            
        self.xtra['loaded'] = scheme+str(domain_index)

def Residue_switch_numbering_scheme(self, scheme, domain_index=0):
    '''
    '''
    # Replace the residue with the numbering scheme with respect to the domain index. The numbering of its own domain is xtra[scheme]
    self.id = ( self.id[0], ) + self.xtra.get(scheme+str(domain_index), self.xtra.get( scheme, self.xtra['pdb'] ) )

def Chain_interface_cys_CA(self, scheme, domain_index=0 ):
    '''
    '''
    self.switch_numbering_scheme( scheme, domain_index )
    try:
        return self[ (' ', inter_cys[self.xtra['domains'][domain_index]][scheme], ' ') ]['CA']
    except KeyError:
        return

def Entity_save(self,filename, select=Select(), write_end=True):
    '''
    '''
    pdbio=PDBIO()
    pdbio.set_structure(self)
    pdbio.save(filename, select=select, write_end=write_end)
        
       
Entity.Entity.switch_numbering_scheme   = Entity_switch_numbering_scheme
Chain.Chain.switch_numbering_scheme     = Chain_switch_numbering_scheme
Chain.Chain.interface_cys_CA            = Chain_interface_cys_CA
Residue.Residue.switch_numbering_scheme = Residue_switch_numbering_scheme
Entity.Entity.save = Entity_save


def compile_remarks( chain, rtype, only_loaded=False ):
    '''
    Compile a remark report that is parsable
    '''
    
    r=[]

    # If only loaded only consider the numbered variable domain on the chain (usually there is only 1)
    if only_loaded:
        r.append( 'CHAIN %s N_%s-V-DOMAINS %d'%( chain.id, rtype , min(1,len(chain.xtra.get('domains',[])) )) )
    else:
        r.append( 'CHAIN %s N_%s-V-DOMAINS %d'%( chain.id, rtype , len(chain.xtra.get('domains',[])) ) )
    
    # Create remark records about the chain types, pairings and germline assignments
    for i in range( len( chain.xtra.get('domains',[]) ) ):
        if only_loaded and i != int(chain.xtra.get('loaded','9')[-1]): continue
        deets = chain.xtra['anarci_details'][i]
        di = 'CHAIN %s DOMAIN %d '%( chain.id, i )
        r.append(di+'TYPE V%s'%chain.xtra['domains'][i]) 
        if chain.xtra['pairing'][i] is not None:
            r.append(di+'PAIREDWITH CHAIN %s DOMAIN %d'%( chain.xtra['pairing'][i][0].id, chain.xtra['pairing'][i][1] ))
        r.append(di+'V_GENE_SPECIES %s'%deets['germlines']['v_gene'][0][0])
        r.append(di+'V_GENE_GERMLINE %s'%deets['germlines']['v_gene'][0][1])
        r.append(di+'V_GENE_GERMLINE_IDENTITY %.2f'%deets['germlines']['v_gene'][1])
        r.append(di+'J_GENE_SPECIES %s'%deets['germlines']['j_gene'][0][0])
        r.append(di+'J_GENE_GERMLINE %s'%deets['germlines']['j_gene'][0][1])
        r.append(di+'J_GENE_GERMLINE_IDENTITY %.2f'%deets['germlines']['j_gene'][1])

    return r        

def compile_seqres( chain ):
    ''' 
    '''
    if 'seqresobj' in chain.xtra:
        return ''.join( [ l%chain.id for l in chain.xtra['seqresobj'].lines ] ).strip()
    return ''

        
def rename_chains( structure ):
    '''
    Rename the first paired receptor chains with (H, L), (B, A). If single domain then only that is renamed.
    Deletes all other chains from each model.
    '''
    for model in structure:
        # Find which chains are the first identified
        keep = [c.id for c,di in model.xtra.get('firstpair',[]) ]
        remove = [ c.id for c in model if c.id not in keep ]

        # Delete unwanted chains
        for ID in remove:
            model.detach_child(ID)

        # Rename wanted chains (detach then add with the new ID)
        lastchain=''
        for chain, di in model.xtra.get('firstpair',[]):
            # Make sure the chain identifier is wrt the first variable domain in the chain. Handles scfvs
            if lastchain == chain.id: continue 

            # Assign the chain identidier with the chain type (Both Kappa and Lambda labelled 'L')
            model.detach_child(chain.id)
            chain.id = chain.xtra['domains'][di].replace('K','L')
            model.add(chain)

            # Update the last chain processed.
            lastchain = chain.id  

def split_scfv( structure ):
    '''
    Split chains with two variable domains into two chains.
    The first has the original identifier. The second has the lower case identifier.
    '''
    for model in structure:
        chains = [c for c in model]
        for chain in chains: 
            if len( chain.xtra.get('domains',[]) ) == 2 and chain.id.lower() not in model: 
                # Biopython builder does not like multiply defined residues (fair enough). Switch back to the original numbering
                # before copying the object as you will have many '0' residues in the first domain.
                scheme = chain.xtra['loaded'][:-1]
                chain.switch_numbering_scheme('pdb')

                # Copy the chain and assign it the lower case identifier
                c_copy = chain.copy()
                c_copy.id = chain.id.lower()
                chain.parent.add(c_copy)                    
        
                # Switch back to the numbering scheme 
                chain.switch_numbering_scheme( scheme, 0 )
                c_copy.switch_numbering_scheme( scheme, 1 )

                # Deal with the pairing changes.
                chain.xtra['pairing'] = [(c_copy,1), None]  
                c_copy.xtra['pairing'] = [None, (chain,0)]
                if chain.id == model.xtra['firstpair'][0][0].id:
                    model.xtra['firstpair'] = [ model.xtra['firstpair'][0], (c_copy,1) ]

                # Reset the region annotations 
                for r in chain: 
                    if r.xtra.get('domain_index','') != 0:
                        r.xtra['region']=''
                for r in c_copy: 
                    if r.xtra.get('domain_index','') != 1:
                        r.xtra['region']=''

# Utility functions
def extract_sequence(chain):
    '''
    Extract the sequence of the chain (amino acids only)
    '''
    strucSeq = ''
    index2ResId = []
    index = 0
    for r in chain:
        a = convert_3_to_1( r.get_resname() )
        if a: # If an amino acid
            index2ResId.append( (index, r.id ) )
            strucSeq += a
            index+=1
    return strucSeq, dict(index2ResId)

def convert_3_to_1(r):
    '''
    Converts amino acid three letter code to one letter codes. '' is returned if a coversion cannot be made.
    '''                     
    if is_aa( r, standard=False ):
        a = protein_letters_3to1.get( r, '' )
        if len( a ) == 3: # Where there are multiple codes used. 
            return self.convert_3_to_1(a)
        return a
    return ''


def pairwise_muscle(seq1, seq2,exact=True):
    """
    Interface with pairwise muscle between two sequences that should be identical.
 
    Try an easy alignment first by checking that one is in the other.
 
    Then if this fails (gaps) use muscle to align the sequences - this should work for seqres and structure sequences with missing atoms.
    
    Use muscle to align.
    
    if exact is True
    The gap open penalty is slightly larger than gap extend to break degeneracy between:
    
    1)garbleabcdefg----
      -a-----bcdefg----
     
    and
     
    2)garbleabcdefg----
      ------abcdefg----
    But they are penalties will still make mismatch very unlikely.
    else
    default muscle
    """

    if not muscle_path:
        raise Exception("Muscle could not be found in the path")    
    
    try_easy  = easy_alignment(seq1, seq2)
    if try_easy:
        return try_easy[0],try_easy[1]
    
    if exact:
        p = subprocess.Popen( [muscle_path,'-gapopen', '-1.001','-gapextend', '-1'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    else:
        p = subprocess.Popen( [muscle_path], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    output = p.communicate('>seq1\n%s\n>seq2\n%s'%(seq1,seq2))
    # check what you have been given.
    result = output[0].split(">")
    # expect it to have 2 entries
    if len(result) == 3:
        seq1_ali = result[1]
        seq2_ali = result[2]
        return "".join( seq1_ali.split("\n")[1:]), "".join( seq2_ali.split("\n")[1:])
    else:
        print("Problem parsing output from muscle: %s"%output[0], file=sys.stderr)

def easy_alignment(seq1, seq2):
    """
    Function to align two sequences by checking if one is in the other.
    This function will conserve gaps.
    """
    assert type(seq1) is str and type(seq2) is str, "Sequences must be strings for easy_alignment" 
    if seq1 in seq2:
        start = seq2.index(seq1)
        seq1_ali = "-"*start + seq1 + "-"*(len(seq2) - start - len(seq1) )
        return seq1_ali, seq2
    elif seq2 in seq1:
        start = seq1.index(seq2)
        seq2_ali = "-"*start + seq2 + "-"*(len(seq1) - start - len(seq2) )
        return seq1, seq2_ali
    else:
        # Can't align them # I return just one value here. 
        return False      


def which(name, flags=os.X_OK):
    """
    Search PATH for executable files with the given name.
   
    On newer versions of MS-Windows, the PATHEXT environment variable will be
    set to the list of file extensions for files considered executable. This
    will normally include things like ".EXE". This fuction will also find files
    with the given name ending with any of these extensions.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
   
    @type name: C{str}
    @param name: The name for which to search.
   
    @type flags: C{int}
    @param flags: Arguments to L{os.access}.
   
    @rtype: C{list}
    @param: A list of the unique full paths to files found, in the
    order in which they were found.
    """
    result = []
    exts = [_f for _f in os.environ.get('PATHEXT', '').split(os.pathsep) if _f]
    path = os.environ.get('PATH', None)
    if path is None:
        return []
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    return uniq(result)

def uniq(seq, idfun=None):
    """
    A function to uniquify a sequence preserving order
    With thanks to http://www.peterbe.com

    @param seq: A sequence to uniquify
    @param idfun: An optional function to use as a key. Like the "key" kwarg in C{sorted}. 
    
    @return: The sequence.
    """
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


def get_alignment_dict(ali1,ali2):
    """
    Get a dictionary which tells you the index in sequence 2 that should align with the index in sequence 1 (key)
    
    ali1:  ----bcde-f---        seq1: bcdef
    ali2:  ---abcd--f---        seq2: abcdf

    alignment_dict={
        0:1,
        1:2,
        2:3,
        4:4
        }
    
    If the index is aligned with a gap do not include in the dictionary.
    e.g  1 in alignment_dict  --> True
    e.g  3 in alignment_dict  --> False
    """
    assert len(ali1)==len(ali2), "aligned sequences must be same lengths (including gaps)"
    alignment_dict={}
    p1=-1
    p2=-1
    for ap in range( len(ali1) ):
        if ali1[ap] != "-" and ali2[ap] != "-":
            assert ali1[ap] == ali2[ap], 'The sequence in the structure did not match that in the provided sequence (seqres entry)'
            p1+=1
            p2+=1
            alignment_dict[p1] = p2
        elif ali1[ap] != "-": 
            p1+=1 
        elif ali2[ap] != "-": 
            p2+=1    
    return alignment_dict


def fix_deletions(chain):
    '''
    Consider the case where there is one or more missing residues flanked by identical residues:
   
    QETGGGGGGNCP
    QETGG-GGGNCP

    The alignment is degenerate. However missing residues can be any one of them

    e.g.
    QETGGGGGGNCP
    QET-GGGGGNCP

    QETGGGGGGNCP
    QETGGGG-GNCP

    This function will check the position of the chain break and fix which positions should be included accordingly

    '''
    # TODO
    pass

def analyse_continuity( chain ):
    '''
    '''
    last = chain.child_list[0]
    i = 1
    for residue in chain.child_list[1:]:
        try:
            if residue.id[0].strip():continue # Don't consider HETATMS. Enforce them to be missed
            if abs(last['C'] - residue['N'] ) > 1.5: # Generous cut-off (average 1.33) for C-N peptide bond.
                return False, i
        except KeyError:
            return False, i
        last = residue
        i+=1
    return True, None

# Check whether muscle is in the path. Alternatively provide this as an input argument from another script. 
if which('muscle'):
    muscle_path = 'muscle'
else:
    muscle_path = ''



# Modified seqres iterator from Biopython
def PdbSeqresIterator(handle):
    """Returns SeqRecord objects for each chain in a PDB file.

    The sequences are derived from the SEQRES lines in the
    PDB file header, not the atoms of the 3D structure.

    Specifically, these PDB records are handled: DBREF, SEQADV, SEQRES, MODRES

    See: http://www.wwpdb.org/documentation/format23/sect3.html
    """
    chains = collections.defaultdict(list)
    metadata = collections.defaultdict(list)
    lines = collections.defaultdict(list)
    for line in handle:
        rec_name = line[0:6].strip()
        if rec_name == 'SEQRES':
            # NB: We only actually need chain ID and the residues here;
            # commented bits are placeholders from the wwPDB spec.
            # Serial number of the SEQRES record for the current chain.
            # Starts at 1 and increments by one each line.
            # Reset to 1 for each chain.
            # ser_num = int(line[8:10])
            # Chain identifier. This may be any single legal character,
            # including a blank which is used if there is only one chain.
            chn_id = line[11]
            # Number of residues in the chain (repeated on every record)
            # num_res = int(line[13:17])
            residues = [seq1(res, custom_map=protein_letters_3to1) for res in line[19:].split()]
            chains[chn_id].extend(residues)
            lines[chn_id].append( line[:11]+'%s'+line[12:] )
        elif rec_name == 'DBREF':
            #  ID code of this entry (PDB ID)
            pdb_id = line[7:11]
            # Chain identifier.
            chn_id = line[12]
            # Initial sequence number of the PDB sequence segment.
            # seq_begin = int(line[14:18])
            # Initial insertion code of the PDB sequence segment.
            # icode_begin = line[18]
            # Ending sequence number of the PDB sequence segment.
            # seq_end = int(line[20:24])
            # Ending insertion code of the PDB sequence segment.
            # icode_end = line[24]
            # Sequence database name.
            database = line[26:32].strip()
            # Sequence database accession code.
            db_acc = line[33:41].strip()
            # Sequence database identification code.
            db_id_code = line[42:54].strip()
            # Initial sequence number of the database seqment.
            # db_seq_begin = int(line[55:60])
            # Insertion code of initial residue of the segment, if PDB is the
            # reference.
            # db_icode_begin = line[60]
            # Ending sequence number of the database segment.
            # db_seq_end = int(line[62:67])
            # Insertion code of the ending residue of the segment, if PDB is the
            # reference.
            # db_icode_end = line[67]
            metadata[chn_id].append({'pdb_id': pdb_id, 'database': database,
                                    'db_acc': db_acc, 'db_id_code': db_id_code})
        # ENH: 'SEQADV' 'MODRES'

    for chn_id, residues in sorted(chains.items()):
        record = SeqRecord(Seq(''.join(residues), generic_protein))
        record.annotations = {"chain": chn_id}
        if chn_id in metadata:
            m = metadata[chn_id][0]
            record.id =  chn_id
            record.name = "%s:%s" % (m['pdb_id'], chn_id)
            record.description = ("%s:%s %s" % (m['database'],
                                                m['db_acc'],
                                                m['db_id_code']))
            record.lines =lines[chn_id]
            for melem in metadata[chn_id]:
                record.dbxrefs.extend([
                    "%s:%s" % (melem['database'], melem['db_acc']),
                    "%s:%s" % (melem['database'], melem['db_id_code'])])
        else:
            record.id = chn_id
            record.lines =lines[chn_id]
        yield record


######## Region annotations #########

regions_in_chothia = {"L1": {"kabat":(24,34),"chothia":(24,34),"contact":(30,36),"north":(24,34)},
                      "L2": {"kabat":(50,56),"chothia":(50,56),"contact":(46,55),"north":(49,56)},
                      "L3": {"kabat":(89,97),"chothia":(89,97),"contact":(89,96),"north":(89,97)},
                      "H1": {"kabat":(31,35),"chothia":(26,32),"contact":(30,35),"north":(23,35)},
                      "H2": {"kabat":(50,65),"chothia":(52,56),"contact":(47,58),"north":(50,58)},
                      "H3": {"kabat":(95,102),"chothia":(95,102),"contact":(93,101),"north":(93,102)} }


_regions = {'imgt':{}}
_regions['imgt']['L'] = _regions['imgt']['H'] = '11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777'
_regions['kabat'] = {}
_regions['kabat']['L']                        = '11111111111111111111111222222222222222223333333333333334444444444444455555555555555555555555555555555555666666666666677777777777'
_regions['kabat']['H']                        = '11111111111111111111111111111111111222223333333333333344444444444444444444555555555555555555555555555555556666666666677777777777'
_regions['chothia'] = {}
_regions['chothia']['L']                      = '11111111111111111111111222222222222222223333333333333334444444444444455555555555555555555555555555555555666666666666677777777777'
_regions['chothia']['H']                      = '11111111111111111111111111222222222223333333333333333333444444445555555555555555555555555555555555555555556666666666677777777777'
_regions['contact'] = {}
_regions['contact']['L']                      = '11111111111111111111111111111111111222222233333333344444444444444444555555555555555555555555555555555555666666666666777777777777'
_regions['contact']['H']                      = '11111111111111111111111111111122222222223333333333344444444444444455555555555555555555555555555555555555666666666666777777777777'
_regions['north'] = {}
_regions['north']['L']                        = '11111111111111111111111222222222222222223333333333333344444444444444455555555555555555555555555555555555666666666666677777777777'
_regions['north']['H']                        = '11111111111111111111111222222222222222223333333333333344444444444455555555555555555555555555555555555555666666666666677777777777'

# For internal use only. These are not direct conversions and are handled heuristically.
_index_to_imgt_state =  {('chothia', 'H'): {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18, 19:
19, 20: 20, 21: 21, 22: 22, 23: 23, 24: 24, 25: 25, 26: 26, 27: 27, 28: 28, 29: 29, 30: 30, 31: 35, 32: 36, 33: 37, 34: 38, 35: 39, 36: 40, 37: 41,
38: 42, 39: 43, 40: 44, 41: 45, 42: 46, 43: 47, 44: 48, 45: 49, 46: 50, 47: 51, 48: 52, 49: 53, 50: 54, 51: 55, 52: 59, 53: 60, 54: 61, 55: 62, 56:
63, 57: 64, 58: 65, 59: 66, 60: 67, 61: 68, 62: 69, 63: 70, 64: 72, 65: 73, 66: 74, 67: 75, 68: 76, 69: 77, 70: 78, 71: 79, 72: 80, 73: 81, 74: 82,
75: 83, 76: 84, 77: 85, 78: 86, 79: 87, 80: 88, 81: 89, 82: 93, 83: 94, 84: 95, 85: 96, 86: 97, 87: 98, 88: 99, 89: 100, 90: 101, 91: 102, 92: 103,
93: 104, 94: 105, 95: 106, 96: 107, 97: 108, 98: 109, 99: 110, 100: 114, 101: 115, 102: 116, 103: 117, 104: 118, 105: 119, 106: 120, 107: 121, 108:
122, 109: 123, 110: 124, 111: 125, 112: 126, 113: 127}, ('kabat', 'H'): {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12,
13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18: 18, 19: 19, 20: 20, 21: 21, 22: 22, 23: 23, 24: 24, 25: 25, 26: 26, 27: 27, 28: 28, 29: 29, 30: 30, 31:
31, 32: 32, 33: 33, 34: 34, 35: 35, 36: 40, 37: 41, 38: 42, 39: 43, 40: 44, 41: 45, 42: 46, 43: 47, 44: 48, 45: 49, 46: 50, 47: 51, 48: 52, 49: 53,
50: 54, 51: 55, 52: 59, 53: 60, 54: 61, 55: 62, 56: 63, 57: 64, 58: 65, 59: 66, 60: 67, 61: 68, 62: 69, 63: 70, 64: 72, 65: 73, 66: 74, 67: 75, 68:
76, 69: 77, 70: 78, 71: 79, 72: 80, 73: 81, 74: 82, 75: 83, 76: 84, 77: 85, 78: 86, 79: 87, 80: 88, 81: 89, 82: 93, 83: 94, 84: 95, 85: 96, 86: 97,
87: 98, 88: 99, 89: 100, 90: 101, 91: 102, 92: 103, 93: 104, 94: 105, 95: 106, 96: 107, 97: 108, 98: 109, 99: 110, 100: 114, 101: 115, 102: 116, 103:
117, 104: 118, 105: 119, 106: 120, 107: 121, 108: 122, 109: 123, 110: 124, 111: 125, 112: 126, 113: 127}, ('imgt', 'H'): {1: 0, 2: 1, 3: 2, 4: 3, 5:
4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17, 19: 18, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23, 25:
24, 26: 25, 27: 26, 28: 27, 29: 28, 30: 29, 31: 30, 32: 31, 33: 32, 34: 33, 35: 34, 36: 35, 37: 36, 38: 37, 39: 38, 40: 39, 41: 40, 42: 41, 43: 42,
44: 43, 45: 44, 46: 45, 47: 46, 48: 47, 49: 48, 50: 49, 51: 50, 52: 51, 53: 52, 54: 53, 55: 54, 56: 55, 57: 56, 58: 57, 59: 58, 60: 59, 61: 60, 62:
61, 63: 62, 64: 63, 65: 64, 66: 65, 67: 66, 68: 67, 69: 68, 70: 69, 71: 70, 72: 71, 73: 72, 74: 73, 75: 74, 76: 75, 77: 76, 78: 77, 79: 78, 80: 79,
81: 80, 82: 81, 83: 82, 84: 83, 85: 84, 86: 85, 87: 86, 88: 87, 89: 88, 90: 89, 91: 90, 92: 91, 93: 92, 94: 93, 95: 94, 96: 95, 97: 96, 98: 97, 99:
98, 100: 99, 101: 100, 102: 101, 103: 102, 104: 103, 105: 104, 106: 105, 107: 106, 108: 107, 109: 108, 110: 109, 111: 110, 112: 111, 113: 112, 114:
113, 115: 114, 116: 115, 117: 116, 118: 117, 119: 118, 120: 119, 121: 120, 122: 121, 123: 122, 124: 123, 125: 124, 126: 125, 127: 126, 128: 127},
('chothia', 'L'): {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17, 19:
18, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23, 25: 24, 26: 25, 27: 26, 28: 27, 29: 28, 30: 35, 31: 36, 32: 37, 33: 38, 34: 39, 35: 40, 36: 41, 37: 42,
38: 43, 39: 44, 40: 45, 41: 46, 42: 47, 43: 48, 44: 49, 45: 50, 46: 51, 47: 52, 48: 53, 49: 54, 50: 55, 51: 56, 52: 57, 53: 65, 54: 66, 55: 67, 56:
68, 57: 69, 58: 70, 59: 72, 60: 73, 61: 74, 62: 75, 63: 76, 64: 77, 65: 78, 66: 81, 67: 82, 68: 83, 69: 84, 70: 85, 71: 86, 72: 87, 73: 88, 74: 89,
75: 90, 76: 91, 77: 92, 78: 93, 79: 94, 80: 95, 81: 96, 82: 97, 83: 98, 84: 99, 85: 100, 86: 101, 87: 102, 88: 103, 89: 104, 90: 105, 91: 106, 92:
107, 93: 108, 94: 109, 95: 114, 96: 115, 97: 116, 98: 117, 99: 118, 100: 119, 101: 120, 102: 121, 103: 122, 104: 123, 105: 124, 106: 125, 107: 126,
108: 127}, ('martin', 'H'): {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 8, 9: 9, 10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15, 16: 16, 17: 17, 18:
18, 19: 19, 20: 20, 21: 21, 22: 22, 23: 23, 24: 24, 25: 25, 26: 26, 27: 27, 28: 28, 29: 29, 30: 30, 31: 35, 32: 36, 33: 37, 34: 38, 35: 39, 36: 40,
37: 41, 38: 42, 39: 43, 40: 44, 41: 45, 42: 46, 43: 47, 44: 48, 45: 49, 46: 50, 47: 51, 48: 52, 49: 53, 50: 54, 51: 55, 52: 59, 53: 60, 54: 61, 55:
62, 56: 63, 57: 64, 58: 65, 59: 66, 60: 67, 61: 68, 62: 69, 63: 70, 64: 72, 65: 73, 66: 74, 67: 75, 68: 76, 69: 77, 70: 78, 71: 79, 72: 83, 73: 84,
74: 85, 75: 86, 76: 87, 77: 88, 78: 89, 79: 90, 80: 91, 81: 92, 82: 93, 83: 94, 84: 95, 85: 96, 86: 97, 87: 98, 88: 99, 89: 100, 90: 101, 91: 102, 92:
103, 93: 104, 94: 105, 95: 106, 96: 107, 97: 108, 98: 109, 99: 110, 100: 114, 101: 115, 102: 116, 103: 117, 104: 118, 105: 119, 106: 120, 107: 121,
108: 122, 109: 123, 110: 124, 111: 125, 112: 126, 113: 127}, ('kabat', 'L'): {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12:
11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17, 19: 18, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23, 25: 24, 26: 25, 27: 32, 28: 33, 29: 34, 30: 35,
31: 36, 32: 37, 33: 38, 34: 39, 35: 40, 36: 41, 37: 42, 38: 43, 39: 44, 40: 45, 41: 46, 42: 47, 43: 48, 44: 49, 45: 50, 46: 51, 47: 52, 48: 53, 49:
54, 50: 55, 51: 56, 52: 57, 53: 65, 54: 66, 55: 67, 56: 68, 57: 69, 58: 70, 59: 72, 60: 73, 61: 74, 62: 75, 63: 76, 64: 77, 65: 78, 66: 81, 67: 82,
68: 83, 69: 84, 70: 85, 71: 86, 72: 87, 73: 88, 74: 89, 75: 90, 76: 91, 77: 92, 78: 93, 79: 94, 80: 95, 81: 96, 82: 97, 83: 98, 84: 99, 85: 100, 86:
101, 87: 102, 88: 103, 89: 104, 90: 105, 91: 106, 92: 107, 93: 108, 94: 109, 95: 114, 96: 115, 97: 116, 98: 117, 99: 118, 100: 119, 101: 120, 102:
121, 103: 122, 104: 123, 105: 124, 106: 125, 107: 126, 108: 127}, ('imgt', 'L'): {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10,
12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17, 19: 18, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23, 25: 24, 26: 25, 27: 26, 28: 27, 29: 28, 30:
29, 31: 30, 32: 31, 33: 32, 34: 33, 35: 34, 36: 35, 37: 36, 38: 37, 39: 38, 40: 39, 41: 40, 42: 41, 43: 42, 44: 43, 45: 44, 46: 45, 47: 46, 48: 47,
49: 48, 50: 49, 51: 50, 52: 51, 53: 52, 54: 53, 55: 54, 56: 55, 57: 56, 58: 57, 59: 58, 60: 59, 61: 60, 62: 61, 63: 62, 64: 63, 65: 64, 66: 65, 67:
66, 68: 67, 69: 68, 70: 69, 71: 70, 72: 71, 73: 72, 74: 73, 75: 74, 76: 75, 77: 76, 78: 77, 79: 78, 80: 79, 81: 80, 82: 81, 83: 82, 84: 83, 85: 84,
86: 85, 87: 86, 88: 87, 89: 88, 90: 89, 91: 90, 92: 91, 93: 92, 94: 93, 95: 94, 96: 95, 97: 96, 98: 97, 99: 98, 100: 99, 101: 100, 102: 101, 103: 102,
104: 103, 105: 104, 106: 105, 107: 106, 108: 107, 109: 108, 110: 109, 111: 110, 112: 111, 113: 112, 114: 113, 115: 114, 116: 115, 117: 116, 118: 117,
119: 118, 120: 119, 121: 120, 122: 121, 123: 122, 124: 123, 125: 124, 126: 125, 127: 126, 128: 127}, ('martin', 'L'): {1: 0, 2: 1, 3: 2, 4: 3, 5: 4,
6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17, 19: 18, 20: 19, 21: 20, 22: 21, 23: 22, 24: 23, 25: 24,
26: 25, 27: 26, 28: 27, 29: 28, 30: 35, 31: 36, 32: 37, 33: 38, 34: 39, 35: 40, 36: 41, 37: 42, 38: 43, 39: 44, 40: 45, 41: 46, 42: 47, 43: 48, 44:
49, 45: 50, 46: 51, 47: 52, 48: 53, 49: 54, 50: 55, 51: 56, 52: 57, 53: 65, 54: 66, 55: 67, 56: 68, 57: 69, 58: 70, 59: 72, 60: 73, 61: 74, 62: 75,
63: 76, 64: 77, 65: 78, 66: 81, 67: 82, 68: 83, 69: 84, 70: 85, 71: 86, 72: 87, 73: 88, 74: 89, 75: 90, 76: 91, 77: 92, 78: 93, 79: 94, 80: 95, 81:
96, 82: 97, 83: 98, 84: 99, 85: 100, 86: 101, 87: 102, 88: 103, 89: 104, 90: 105, 91: 106, 92: 107, 93: 108, 94: 109, 95: 114, 96: 115, 97: 116, 98:
117, 99: 118, 100: 119, 101: 120, 102: 121, 103: 122, 104: 123, 105: 124, 106: 125, 107: 126, 108: 127}}

wolfguy_indexdiv50_to_region = {'H': [ 'fwh1', 'cdrh1','fwh2', 'cdrh2','fwh3', 'cdrh3', 'fwh4'], 
                                'L': [ 'fwl1', 'cdrl1','fwl2', 'cdrl2','fwl3', 'cdrl3', 'fwl4'] }

_reg_one2three = { "1":"fw%s1","2":"cdr%s1","3":"fw%s2","4":"cdr%s2","5":"fw%s3","6":"cdr%s3","7":"fw%s4" }

def get_region( position, chain, numbering_scheme="chothia", definition="chothia" ):
    """
    Get the region in which the position belongs given the chain, numbering scheme and definition.

    **Note** this function does not know about insertions on the sequence. Therefore, it will get the region annotation
    wrong when using non-equivalent scheme-definitions. 

    To get around this please use the annotate_regions function which implements heuristics to get the definition correct
    in the scheme.

    """
    index, insertion = position
    chain=chain.upper()

    # Horrible exception cases revolving around the kabat scheme/definition and cdr h1
    if definition == "kabat":  
        if numbering_scheme == "kabat" and chain == "H" and 31 <= index < 36: # Kabat scheme kabat definition.
            if index == 35:
                if insertion in " AB": # Position 31 to 35B
                    return "cdrh1"
                else:
                    return "fwh2" # 31C would be framework.                   
            else:
                return "cdrh1"
    if numbering_scheme == "kabat": # Kabat numbering, chothia or imgt definitions.
        if definition=="chothia" and chain == "H" and 33 <= index < 36:
            return "fwh2"
        elif definition=="imgt" and chain == "H" and 34 <= index < 36:
            return "fwh2"

    if numbering_scheme == "wolfguy" or definition == "wolfguy":
        assert definition == "wolfguy" and numbering_scheme == "wolfguy", "The wolfguy numbering scheme must be used with the wolfguy CDR definition"
        if chain == 'H':
            if index > 411: return ""
            r = int(index/50)-2
        elif chain == 'L':
            if index > 810: return ""
            r = int(index/50)-10
        try:
            return wolfguy_indexdiv50_to_region[chain][r]
        except IndexError:
            return ""

    try:
        return _reg_one2three[_regions[definition][chain][  _index_to_imgt_state[ (numbering_scheme, chain) ][index] ]]%chain.lower()
    except KeyError:
        return "?"

# Heuristics to locate the different definitions in different numbering scheme
# TODO
# All definitions in the Aho scheme
# All definitions in wolfguy scheme (apart from wolfguy)
def annotate_regions(numbered_sequence, chain,numbering_scheme="chothia",definition="chothia"):
    """
    Given a numbered sequence annotate which region each residue belongs to.
    
    The numbering scheme can be one chothia, kabat, imgt or martin
    The definition can be chothia, kabat, imgt, north or contact.
    
    Contact definition cannot be used with the kabat numbering scheme.

    If possible, use the corresponding numbering scheme and definition.

    This function automates the heuristics recognise different definitions in each scheme. However, 
    some of the conversions are non-trivial. 

    """
    # In some cases there is not a direct equivalence between numbering schemes. Therefore, a CDR definition
    # may not be easily translatable from scheme to scheme. This is the case for:
    #   - H1 Kabat   definition in the IMGT    scheme   # 1
    #   - H1 Chothia definition in the Kabat   scheme   # 2
    #   - H1 IMGT    definition in the Kabat   scheme   # 3
    #   - H1 Chothia definition in the Kabat   scheme   # 4
    #   - L2 IMGT    definition in the Kabat   scheme   # 4
    #   - L2 IMGT    definition in the Chothia scheme   # 5
    #   - L2 IMGT    definition in the Martin  scheme   # 6

    # Below are heuristics to allow the conversion. They have been developed so that the CDR sequence extracted
    # in all schemes is the same as the CDR sequence in the scheme the definition was originally defined. e.g.
    # Chothia for Chothia, Kabat for Kabat, IMGT for IMGT. Note that the Contact definition *cannot* be defined
    # in Kabat numbering. 

    # Find which additional positions should be considered in the CDR, given the sequence.
    assert numbering_scheme in ['chothia','martin','kabat','imgt'], 'Unimplemented'
    additional_positions = {}
    excluded_positions = {}
    c = chain.lower()

    numdict = dict( numbered_sequence )

    cdr_acceptors = { 1:Accept(numbering_scheme=numbering_scheme, definition=definition),
                      2:Accept(numbering_scheme=numbering_scheme, definition=definition),
                      3:Accept(numbering_scheme=numbering_scheme, definition=definition) }

    cdr_acceptors[1].set_regions( ["cdr%s1"%c] )
    cdr_acceptors[2].set_regions( ["cdr%s2"%c] )
    cdr_acceptors[3].set_regions( ["cdr%s3"%c] )

    # Heavy chain 
    if chain == "H":
        if numbering_scheme == "imgt" and definition == "kabat": # IMGT scheme / Kabat definition
            # 1
            # Count the positions 31-34 inclusive. These become insertions on 35 in Kabat.
            ins = 0 
            for i in range( 31, 35 ): 
                if (i, " ") in numdict: 
                    if numdict[ (i, " ") ] != "-": ins +=1
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (33, a) in numdict: 
                    if numdict[ (33, a) ] != "-": ins +=1
                else: break
            # If insertions would occur on 35 in the Kabat scheme, extend the IMGT numbering back until you hit 31. Then put on the insertions at 32
            if ins: # Add postions backwards based on the number of insertions starting with H35.
                    # IMGT insertions go on 33 
                cdr_acceptors[1].add_positions( ([ (35, " "), (34, " "), (33, " "), (32, " ") ]+[(33, a) for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"])[:ins], "H" )
                # If there more than two insertions we have to start removing positions from the definition. See # 4 below for more of an explanation
                if ins >2:
                    cdr_acceptors[1].exclude_positions( ([(40," "),(39," "),(38," "),(37," "),(36," "),(35," "),(34," ")]+[(33, a) for a in "FGHIJKLMNOPQRSTUVWXYZ"])[:ins-2], "H" )
        elif numbering_scheme == "kabat" and definition in ["chothia","imgt"]: # Kabat scheme / Chothia or IMGT definition
            # 2, 3
            # Count the insertions on 35. These happen outside the CDR definitions of IMGT or Chothia.
            # They therefore are missed by a straight conversion. The equivalence is dependent on the 
            # number of insertions at this position.

            # Count the number of insertions on H35.
            ins = 0
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (35, a) in numdict: 
                    if numdict[(35,a)] != "-":
                        ins +=1
                else: break
            if ins: # Add postions based on the number of insertions starting with H33 for the Chothia definition. H34 for the IMGT definition.
                if definition == "chothia":
                    cdr_acceptors[1].add_positions( ([(33, " "),(34," ")]+[(35, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"])[:ins], "H" )
                elif definition == "imgt":
                    cdr_acceptors[1].add_positions( ([(34, " ")]+[(35, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"])[:ins] , "H" )
        elif numbering_scheme in ["chothia", "martin"] and definition == "kabat":
            # 4
            # Count the insertions on 31. If there are more than two, exclude back from 35. 
            # e.g. if we have the below. Kabat and Chothia numbered Kabat CDRs
            # [((31, ' '), 'P'), ((32, ' '), 'A'), ((33, ' '), 'P'), ((34, ' '), 'E'), ((35, ' '), 'H'), ((35, 'A'), 'F'), ((35, 'B'), 'I')]
            # [((31, ' '), 'P'), ((31, 'A'), 'A'), ((31, 'B'), 'P'), ((31, 'C'), 'E'), ((32, ' '), 'H'), ((33, ' '), 'F'), ((34, ' '), 'I')] 
            # Kabat definition will have a max of length 7 as insertions are on 35 and the definition ends at 35B.
            # In Chothia numbering the insertions are at 31. Hence when more than 2 insertions occur Chothia 35 falls off the end. Chothia 34
            # becomes the end of the Kabat CDR. Similar thing happens with IMGT. 
            # Count the number of insertions on H31.
            ins = 0
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (31, a) in numdict: 
                    if numdict[(31,a)] != "-":
                        ins +=1
                else: break
            if ins > 2:
                cdr_acceptors[1].exclude_positions(  ([ (35," "),(34," "),(33," "), (32," ")]+[(31, a) for a in "GHIJKLMNOPQRSTUVWXYZ"])[:ins-2], "H" )

    # Light chain
    if chain == "L":
        # 5,6,7
        if numbering_scheme in ["kabat", "chothia", "martin"] and definition == "imgt": # Kabat-like schemes and IMGT definition
            # Count the number of insertions on L54.
            ins = 0
            for a in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
                if (54, a) in numdict: 
                    if numdict[(54,a)] != "-":
                        ins +=1
                else: break
            if ins: # Add positions based on the number of insertions starting with L53. 
                    # IMGT definition with no insertions ends and 52 in Kabat/Chothia/Martin.
                extensions = [ (53," ") ] + [ (54, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
                cdr_acceptors[2].add_positions( ([ (53," ") ] + [ (54, a) for a in " ABCDEFGHIJKLMNOPQRSTUVWXYZ"])[:ins], "L")

    
    fw_regions = [ "fw%s1"%c,"fw%s2"%c,"fw%s3"%c,"fw%s4"%c ]
    fw_region = "fw%s1"%c
    region_annotations = []
    cterm = max( _index_to_imgt_state[ (numbering_scheme, chain ) ].keys() )
    for r, a in numbered_sequence:
        if cdr_acceptors[1].accept( r, chain ):
            region_annotations.append( (r, a, "cdr%s1"%c))
            fw_region = "fw%s2"%c
        elif cdr_acceptors[2].accept( r, chain ):
            region_annotations.append( (r, a, "cdr%s2"%c))
            fw_region = "fw%s3"%c
        elif cdr_acceptors[3].accept( r, chain ):
            region_annotations.append( (r, a, "cdr%s3"%c))
            fw_region = "fw%s4"%c
        elif r[0] <= cterm and r[0] > 0: # Anything out of the variable region is not assigned a region i.e. ''
            region_annotations.append( (r, a, fw_region ) )
        else:
            region_annotations.append( (r, a, '' ) )

    return region_annotations


class Accept:
    """
    A class to select which positions should be compared.
    """
    _defined_regions = ["fwh1", "fwh2", "fwh3", "fwh4", "fwl1", "fwl2", "fwl3", "fwl4", "cdrh1", "cdrh2", "cdrh3", "cdrl1", "cdrl2", "cdrl3"]
    _macro_regions= { "hframework":set(["fwh1", "fwh2", "fwh3", "fwh4"]),
                      "hcdrs"     :set(["cdrh1", "cdrh2", "cdrh3"]),
                      "lframework":set(["fwl1", "fwl2", "fwl3", "fwl4"]),    
                      "lcdrs"     :set(["cdrl1", "cdrl2", "cdrl3"]) }
    _macro_regions.update( { "framework": _macro_regions["hframework"] | _macro_regions["lframework"],
                             "cdrs": _macro_regions["hcdrs"] | _macro_regions["lcdrs"],
                             "vh": _macro_regions["hcdrs"] | _macro_regions["hframework"],
                             "vl": _macro_regions["lcdrs"] | _macro_regions["lframework"],
                            })
                           
    _macro_regions.update( { "fv": _macro_regions["vh"] | _macro_regions["vl"] } )

    _macro_positions = {}

    def __init__(self, numbering_scheme="chothia", definition="chothia", NOT=False):
        self.NOT = NOT
        self.set_regions()
        self.positions={"H":set(),"L":set()}
        self.numbering_scheme = numbering_scheme
        self.definition = definition
        self.exclude={"H":set(),"L":set()}

    def set_regions( self, regions=[] ):
        """
        Set the regions to be used. Will clear anything added using add regions.
        """
        if self.NOT:
            self.regions= self._macro_regions[ "fv" ]
        else:   
            self.regions = set()
        self.add_regions( regions )

    def add_regions(self, regions):
        """
        Add regions to the selection. 
        """
        for region in regions:
            region = region.lower()
            if region in self._defined_regions:
                if self.NOT:
                    self.regions = self.regions - set([ region ])
                else:
                    self.regions.add( region )
            elif region in self._macro_regions:
                if self.NOT:
                    self.regions = self.regions - self._macro_regions[region]
                else:
                    self.regions = self.regions | self._macro_regions[region]
            elif region in self._macro_positions: # e.g. interface positions
                raise AssertionError
            else:
                raise AssertionError

    def add_positions( self, positions, chain ):
        for position in positions:
            index, insertion = position
            self.positions[chain].add( (index, insertion) )

    def exclude_positions( self, positions, chain ):
        for position in positions:
            index, insertion = position
            self.exclude[chain].add( (index, insertion) )

    def accept(self, position, chain):
        if position in self.exclude[chain]: return
        if get_region(position, chain, self.numbering_scheme, self.definition) in self.regions or position in self.positions[chain]:
            return 1

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(prog="ImmunoPDB", description=description, epilog=epilogue,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument( '-i', type=str, help="A structure to be numbered", dest="inputstructure")
    parser.add_argument( '-o', type=str, default=False, help="The output file to use. Default is stdout", dest="outfile")
    parser.add_argument( '--scheme','-s', type=str, choices=scheme_names+['pdb'], default="imgt", help="Which numbering scheme should be used. i, k, c, m, w and a are shorthand for IMGT, Kabat, Chothia, Martin (Extended Chothia), Wolfguy and Aho respectively. Default IMGT. Use pdb to retain the numbering but get the annotations as remarks", dest="scheme")
    parser.add_argument( '--receptor','-r', type=str, choices=["ig","tr"], default='ig', help="Choose whether to number Antibody (ig) or TCR (tr) domains.", dest="receptor")
    parser.add_argument( '--rename', action='store_true',default=False,help='Rename the receptor chains with H and L (ig) or B and A (tr). Only receptor chains output. First pair identified used.',dest='rename')
    parser.add_argument( '--fvonly', action='store_true',default=False,help='Only output Fv regions.',dest='fvonly')
    parser.add_argument( '--splitscfv', action='store_true',default=False,help='When they are found split single chain fvs into two seperate chains (fvonly becomes true)',dest='splitscfv')
    parser.add_argument( '--warnings', action='store_true',default=False,help='Report warnings about missing residues',dest='warnings')        
    args = parser.parse_args()
    
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)


    # Work arounds for either returning pdb scheme or where I have not yet implemented the regions in the scheme.
    # Assign the regions in imgt then renumber back in the chosen scheme
    # This is slow for the second case (double numbering).
    switchback=False    
    if args.scheme == 'pdb' or ((args.fvonly or args.splitscfv) and args.scheme in ['aho','wolfguy']): 
        switchback = args.scheme
        args.scheme = 'imgt'

    if args.scheme != 'pdb':
        args.scheme = scheme_short_to_long[args.scheme.lower()]

    if args.receptor == 'ig':
        sparser = AntibodyPDBParser(QUIET=True, scheme=args.scheme, warnings=args.warnings)
    elif args.receptor == 'tr':
        if args.scheme not in ['imgt','i','aho','a']: 
            print('Only imgt or aho schemes can be applied to tcrs', file=sys.stderr)
            sys.exit(1)
        sparser = TcrPDBParser(QUIET=True, scheme=args.scheme, warnings=args.warnings)

    name = os.path.splitext(os.path.split(args.inputstructure)[1])[0]
    try:        
        structure = sparser.get_structure(name, args.inputstructure)
    except IOError:
        print('File %s could not be opened'%args.inputstructure, file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print('%s could not be parsed:'%args.inputstructure, e, file=sys.stderr)
        sys.exit(1)

    # Switch back to the numbering scheme of choice. 
    if switchback:
        structure.switch_numbering_scheme(switchback)
        args.scheme=switchback

    # Output the structure according to the input options. 
    try:
        if args.outfile:
            of = open(args.outfile,'w')
        else:
            of = sys.stdout
        print('REMARK   6 SCHEME %s'%args.scheme.upper(), file=of)
        print('REMARK   6 ANARCI TYPE, PAIRING AND ASSIGNED GERMLINE DETAILS', file=of)
        if args.splitscfv: split_scfv( structure )         
        if args.rename: rename_chains( structure )
        for chain in structure.get_chains():
            print('REMARK   6 '+'\nREMARK   6 '.join(compile_remarks(chain, args.receptor.upper(), only_loaded=args.splitscfv)), file=of)
        for chain in structure.get_chains():
            sr = compile_seqres( chain )
            if sr:
                print(sr, file=of)
        if args.fvonly:
            structure.save(of,select=SelectFv())
        elif args.splitscfv:
            structure.save(of,select=SelectFvScFv())
        else:
            structure.save(of)
        if args.outfile:
            of.close()
    except IOError:
        print('Could not write to file %s'%(args.outfile), file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit(1)
    sys.exit(0)

