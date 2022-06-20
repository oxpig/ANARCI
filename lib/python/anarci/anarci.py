#    ANARCI - Antibody Numbering and Antigen Receptor ClassIfication
#    Copyright (C) 2016 Oxford Protein Informatics Group (OPIG)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.#
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
ANARCI - Antigen Receptor Numbering And ClassIfication

Oxford Protein Informatics Group (OPIG). 2015-17

ANARCI performs alignments of sequences to databases of Hidden Markov Models (HMMs).
Those that align with a significant score are classified by species and chain type.
They are then numbered with a scheme of the user's choosing. 

Currently implemented schemes: 
    IMGT
    Chothia (IGs only)
    Kabat (IGs only)
    Martin / Enhanced Chothia (IGs only)
    AHo 
    Wolfguy (IGs only)

Currently recognisable species (chains):
    Human (heavy, kappa, lambda, alpha, beta)
    Mouse (heavy, kappa, lambda, alpha, beta)
    Rat (heavy, kappa, lambda)
    Rabbit (heavy, kappa, lambda)
    Pig (heavy, kappa, lambda)
    Rhesus Monkey (heavy, kappa)
    
Notes:
 o Use assign_germline to get a better species assignment
 o Each scheme has been implemented to follow the published specification as closely as possible. However, in places some schemes
   do not specifiy where insertions should be placed (e.g. imgt FW3). In these cases the HMM alignment is used. This can give rise
   to inserted positions that were not described by the respective paper. 
 o AHo is implemented heuristically based on chain type. If one grafted a foreign CDR1 loop onto, say, a VH domain, it will be 
   numbered as if it is a CDRH1 loop. 
    

'''

import os
import sys
import gzip
import math
import traceback
from functools import partial
from textwrap import wrap
from itertools import groupby, islice
from multiprocessing import Pool

import pyhmmer

# Import from the schemes submodule
from .schemes import *
from .germlines import all_germlines
    
all_species = list(all_germlines['V']['H'].keys())

amino_acids = sorted(list("QWERTYIPASDFGHKLCVNM"))
set_amino_acids = set(amino_acids)
anarci_path  = os.path.split(__file__)[0]

scheme_short_to_long = { "m":"martin", "c":"chothia", "k":"kabat","imgt":"imgt", "kabat":"kabat", "chothia":"chothia", "martin":"martin", "i":"imgt", "a":"aho","aho":"aho","wolfguy":"wolfguy", "w":"wolfguy"}

scheme_names = list(scheme_short_to_long.keys()) 
chain_type_to_class = {"H":"H", "K":"L", "L":"L", "A":"A", "B":"B", "G":"G", "D":"D"}

HMM_path =  os.path.join( anarci_path, "dat", "HMMs" )

all_reference_states = list(range( 1, 129)) # These are the IMGT reference states (matches)

class HMMscanError(Exception):
    def __init__(self, message):
        # Call the base class constructor with the parameters it needs
        super(HMMscanError, self).__init__(message)

## Utility functions ##
def read_fasta(filename):
    """
    Read a sequence file and parse as description, string 
    """
    return [ r for r in fasta_iter(filename) ]

def fasta_iter(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    https://www.biostars.org/p/710/
    """
    if fasta_name.endswith( '.gz' ): # IOError raised upon iteration if not a real gzip file.
        fh = gzip.open(fasta_name)
    else:
        fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        #header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def write_fasta(sequences, f):
    """
    Write a list of sequences to file. 
    
    should be a list of name, sequence tuples
    
    f should be an open file
    """
    for name, sequence in sequences:
        print(">%s"%name, file=f)
        print('\n'.join(['\n'.join(wrap(block, width=80)) for block in sequence.splitlines()]), file=f)
    
    
def validate_sequence(sequence):
    """
    Check whether a sequence is a protein sequence or if someone has submitted something nasty.
    """
    assert len(sequence) < 10000, "Sequence too long."
    assert not (set( sequence.upper() ) - set_amino_acids ), "Unknown amino acid letter found in sequence: %s"% ", ".join(list((set( sequence.upper() ) - set_amino_acids )))    
    return True

def validate_numbering(xxx_todo_changeme, name_seq=[]):
    """
    Wrapper to do some basic validation of the numbering.
    
    Further validation could be done but at the moment we just check that the numbering indices are incremental (they should be)
    """
    (numbering, start, end) = xxx_todo_changeme
    name, seq = name_seq
    last = -1
    nseq=""

    for (index, _), a in numbering:
        assert index >= last, "Numbering was found to decrease along the sequence %s. Please report."%name
        last = index
        nseq += a.replace("-","")

    assert nseq in seq.replace("-",""), "The algorithm did not number a contiguous segment for sequence %s. Please report"%name

    return numbering, start, end

def grouper(n, iterable):
    '''
    Group entries of an iterable by n
    '''
    it = iter(iterable)
    def take():
        while 1:
            yield list( islice(it,n) )
    return iter(take().__next__, [] )

def anarci_output(numbered, sequences, alignment_details, outfile, sequence_id=None, domain_id=None):
    """
    Outputs to open file

    If sequence_id is specified as an integer then only this sequence will be printed. 
    Otherwise all sequences will be printed.

    If domain_id is specified as an integer then only this domain will be printed. 
    Otherwise all domains will be printed.

    If domain_id is specified then sequence_id must also be specified. 
    """       
    assert (sequence_id is not None) or (sequence_id is None and domain_id is None), "If domain_id is specified, sequence_id must also be specified."
    for i in range(len(numbered)):
        if sequence_id is None:
            print("# %s"%sequences[i][0], file=outfile) # print the name
        if numbered[i] is not None:
            if sequence_id is not None:
                if i != sequence_id: continue
            print("# ANARCI numbered", file=outfile)
            for j in range( len(numbered[i])): # Iterate over domains
                if domain_id is not None:
                    if j != domain_id: continue
                print("# Domain %d of %d"%(j+1, len(numbered[i]) ), file=outfile)
                print("# Most significant HMM hit", file=outfile)
                print("#|species|chain_type|e-value|score|seqstart_index|seqend_index|", file=outfile)
                alignment_details[i][j]["evalue"] = str( alignment_details[i][j]["evalue"] )
                print("#|%s|%s|%s|%.1f|%d|%d|"%tuple( [alignment_details[i][j][field] for field in 
                                                                     ["species","chain_type","evalue","bitscore"]] 
                                                                   +[ numbered[i][j][1], numbered[i][j][2]] ), file=outfile)
                
                if 'germlines' in alignment_details[i][j]:
                    print('# Most sequence-identical germlines', file=outfile)
                    print('#|species|v_gene|v_identity|j_gene|j_identity|', file=outfile)
                    (species, vgene), vid =alignment_details[i][j]['germlines'].get('v_gene', [['','unknown'],0])
                    if vgene is None:
                        vgene, vid = 'unknown', 0
                    (_,jgene), jid =alignment_details[i][j]['germlines'].get('j_gene', [['','unknown'],0])
                    if jgene is None:
                        jgene, jid = 'unknown', 0
                    print('#|%s|%s|%.2f|%s|%.2f|'%(species, vgene, vid, jgene, jid ), file=outfile)	
                chain_type = chain_type_to_class[  alignment_details[i][j]["chain_type"] ]
                print("# Scheme = %s"%alignment_details[i][j]["scheme"], file=outfile)
                if len( numbered[i][j][0] ) == 0:
                    print("# Warning: %s scheme could not be applied to this sequence."%alignment_details[i][j]["scheme"], file=outfile)            
                for (index, insertion), aa in numbered[i][j][0]:
                    print(chain_type, ("%d"%index).ljust(5), insertion, aa, file=outfile)
        print("//", file=outfile)

def csv_output(sequences, numbered, details, outfileroot):
    '''
    Write numbered sequences to csv files. A csv file is written for each chain type.

    Kappa and Lambda chains are written to the same file

    The sequences will written aligned to the numbering scheme. Gaps in the sequences with respect to the alignment are written
    as a '-'

    @param sequences: List of name, sequence tuples    
    @param numbered: Numbered sequences in the same order as the sequences list. 
    @param details: List of alignment details in the same order as the sequences list.
    @param outfileroot: The file path for csv files to write. _<chain_type>.csv will be appended to this.
    '''
    
    chain_types = {}        
    pos_ranks = {}
    all_pos = {}
    _lc = {'K':'KL','L':'KL'}


    # Divide the set into chain types and find how to order the numbering for each type.
    for i in range( len(sequences) ): # Iterate over entries
        if numbered[i] is None: continue

        for j in range(len(numbered[i])): # Iterate over domains.
            # Record the chain type index
            c = details[i][j]['chain_type']
            c = _lc.get(c, c) # Consider lambda and kappa together.
            chain_types.setdefault( c, [] ).append( (i,j) ) 
            if c not in pos_ranks:
                pos_ranks[c] = {}
                all_pos[c] = set()

            # Update the insertion order for the scheme. i.e. is it A B C or C B A (e.g. imgt 111 and 112 repectively)
            l = -1 
            r = 0 
            for p, _ in numbered[i][j][0]:
                if p[0] != l:
                    l = p[0]
                    r = 0
                else:
                    r +=1
                pos_ranks[c][p] = max( r, pos_ranks[c].get( p, r ) )
                all_pos[c].add( p )

    # Write a new file for each chain type. Kappa and lambda are written together as light chains. 
    for cts in ['H','KL','A','B','G','D']:
        if cts in chain_types:
            with open( outfileroot + '_%s.csv'%cts, 'w' ) as out:

                # Sort the positions by index and insertion order
                positions = sorted( all_pos[cts], key = lambda p: (p[0], pos_ranks[cts][p]) )

                # Header line
                fields = ['Id','domain_no','hmm_species','chain_type','e-value','score','seqstart_index','seqend_index',
                          'identity_species','v_gene','v_identity','j_gene','j_identity']
                fields += [ ('%d%s'%(p)).strip() for p in positions ]
                print(','.join( fields ), file=out)

                # Iterate over the domains identified
                for i,j in chain_types[cts]:
                    line = [ sequences[i][0].replace(',',' '),      
                             str(j),
                             details[i][j].get('species',''),
                             details[i][j].get('chain_type',''),
                             str(details[i][j].get('evalue','')),
                             str(details[i][j].get('bitscore','')),
                             str(numbered[i][j][1]),
                             str(numbered[i][j][2]),
                             details[i][j].get('germlines',{}).get( 'v_gene',[['',''],0] )[0][0],
                             details[i][j].get('germlines',{}).get( 'v_gene',[['',''],0] )[0][1],
                             '%.2f'%details[i][j].get('germlines',{}).get( 'v_gene',[['',''],0] )[1],
                             details[i][j].get('germlines',{}).get( 'j_gene',[['',''],0] )[0][1],
                             '%.2f'%details[i][j].get('germlines',{}).get( 'j_gene',[['',''],0] )[1] ]

                    # Hash the numbering. Insertion order has been preserved in the positions sort.
                    d = dict( numbered[i][j][0] )
                    line += [ d.get(p,'-') for p in positions ]

                    assert len( line ) == len( fields )
                    print(','.join( line ), file=out)



## Parsing and recognising domain hits from hmmscan ##
def _domains_are_same(dom1, dom2):
    """
    Check to see if the domains are overlapping.
    @param dom1: 
    @param dom2: 

    @return: True or False  
    """
    dom1, dom2 = sorted( [dom1, dom2], key=lambda dom: dom.alignment.target_from  )
    if dom2.alignment.target_from-1 >= dom1.alignment.target_to:
        return False
    return True


def _parse_hmmer_query(hsps, seq_len, hmmer_species=None):
    """
    
    @param hsps: List of Domain objects from pyhmmer (corresponding to BioPython hsps)
    @param seq_len: sequence length

    The function will identify multiple domains if they have been found and provide the details for the best alignment for each domain.
    This allows the ability to identify single chain fvs and engineered antibody sequences as well as the capability in the future for identifying constant domains. 

    """
    hit_table = [ ['id', 'description', 'evalue', 'bitscore', 'bias', 
                    'query_start', 'query_end' ] ]

    # Find the best hit for each domain in the sequence.

    top_descriptions, domains,state_vectors = [], [], []
    
    if hsps: # We have some hits
        # If we have specified a species, check to see we have hits for that species
        # Otherwise revert back to using any species
        if hmmer_species:
            hit_correct_species = []
            for hsp in hsps:
                for species in hmmer_species:
                    hit_id = hsp.alignment.hmm_name.decode()
                    if hit_id.startswith(species):
                        hit_correct_species.append(hsp)

            if hit_correct_species:
                hsp_list = hit_correct_species
            else:
                print("Limiting hmmer search to species %s was requested but hits did not achieve a high enough bitscore. Reverting to using any species" %(hmmer_species))
                hsp_list = hsps
        else:
            hsp_list = hsps

        for hsp in sorted(hsp_list, key=lambda x: x.i_evalue): # Iterate over the matches of the domains in order of their e-value (most significant first)
            new=True
            for i in range( len(domains) ): # Check to see if we already have seen the domain
                if _domains_are_same( domains[i], hsp ):
                    new = False
                    break
            ali = hsp.alignment
            hit_table.append( [ ali.hmm_name.decode(), ali.hmm_accession.decode(), hsp.i_evalue, hsp.score, hsp.bias, ali.target_from-1, ali.target_to] )
            if new: # It is a new domain and this is the best hit. Add it for further processing.
                domains.append( hsp )
                top_descriptions.append(  dict( list(zip(hit_table[0], hit_table[-1])) ) ) # Add the last added to the descriptions list.

        # Reorder the domains according to the order they appear in the sequence.         
        ordering = sorted( list(range(len(domains))), key=lambda x: domains[x].alignment.target_from)
        domains = [ domains[_] for _ in ordering ]
        top_descriptions = [ top_descriptions[_] for _ in ordering ]         
   
    ndomains = len( domains )
    for i in range(ndomains): # If any significant hits were identified parse and align them to the reference state.
        species, chain = top_descriptions[i]["id"].split("_")
        state_vectors.append( _hmm_alignment_to_states(domains[i], ndomains, seq_len, i) ) # Alignment to the reference states.
        top_descriptions[i][ "species"] = species # Reparse
        top_descriptions[i][ "chain_type"] = chain
        top_descriptions[i][ "query_start"] = state_vectors[-1][0][-1] # Make sure the query_start agree if it was changed

    return hit_table, state_vectors, top_descriptions


def _hmm_alignment_to_states(hsp, n, seq_length, order):
    """
    Take a hit hsp and turn the alignment into a state vector with sequence indices
    """
    ali = hsp.alignment

    # Extract the start an end points of the hmm states and the sequence
    # These are python indices i.e list[ start:end ] and therefore start will be one less than in the text file
    _hmm_start = ali.hmm_from - 1
    _hmm_end = ali.hmm_to
     
    _seq_start = ali.target_from - 1
    _seq_end = ali.target_to

    # Extract the full length of the HMM hit
    species, ctype = ali.hmm_name.decode().split('_')
    _hmm_length = get_hmm_length( species, ctype )

    hmm_seq = ali.hmm_sequence.upper()
    target_seq = ali.target_sequence

    # Handle cases where there are n terminal modifications.
    # In most cases the user is going to want these included in the numbered domain even though they are not 'antibody like' and 
    # not matched to the germline. Only allow up to a maximum of 5 unmatched states at the start of the domain
    # Adds a bug here if there is a very short linker between a scfv domains with a modified n-term second domain
    # Thus this is only done for the first identified domain ( hence order attribute on hsp )
    if order == 0 and _hmm_start and _hmm_start < 5:
        n_extend = _hmm_start 
        if _hmm_start > _seq_start:
            n_extend = min( _seq_start , _hmm_start - _seq_start )
        hmm_seq = 'X'*n_extend + hmm_seq
        target_seq = 'X'*n_extend + target_seq
        _seq_start = _seq_start - n_extend
        _hmm_start = _hmm_start - n_extend

    # Handle cases where the alignment should be extended to the end of the j-element
    # This occurs when there a c-terminal modifications of the variable domain that are significantly different to germline
    # Extension is only made when half of framework 4 has been recognised and there is only one domain recognised.
    if n==1 and _seq_end < seq_length and (123 < _hmm_end < _hmm_length): # Extend forwards
        n_extend = min( _hmm_length - _hmm_end, seq_length - _seq_end )
        hmm_seq = hmm_seq + 'X'*n_extend
        target_seq = target_seq + 'X'*n_extend
        _seq_end = _seq_end + n_extend
        _hmm_end = _hmm_end + n_extend

    assert len(hmm_seq) == len(target_seq), \
        'Unexpected mismatch in aligned HMM sequence, {} != {}'.format(len(hmm_seq), len(target_seq))

    # Generate state vector from aligned HMM sequence and input sequence
    # Each item in the list corresponds to: (hmm_state, state_type), sequence_index
    # hmm_state starts from 1, sequence_index starts from 0
    # For example: (111, 'm'), 110
    state_vector = []
    hmm_state = _hmm_start + 1 # IMGT HMM states go from 1 to 128
    sequence_index = _seq_start
    for hmm_aa, target_aa in zip(hmm_seq, target_seq):
        state_type = 'i' if hmm_aa == '.' else ('d' if target_aa == '-' else 'm')
        state_vector.append(((hmm_state, state_type), sequence_index if state_type != 'd' else None))
        if state_type == 'm' or state_type == 'i':
            sequence_index += 1
        if state_type == 'm' or state_type == 'd':
            hmm_state += 1

    assert hmm_state == _hmm_end + 1, \
        'Unexpected mismatch in HMM end {} != {}'.format(hmm_state, _hmm_end + 1)
    assert sequence_index == _seq_end, \
        'Unexpected mismatch in input end {} != {}'.format(sequence_index, _seq_end)

    return state_vector


def run_pyhmmer(sequence_list, hmm_database="ALL", ncpu=None, bit_score_threshold=80, hmmer_species=None):
    """
    Run the sequences in sequence list against a precompiled hmm_database.

    Those sequence that have a significant hit with a bit score over a threshold will
    be recognised and an alignment given. The alignment will be used to number the
    sequence.

    @param sequence_list: a list of (name, sequence) tuples. Both are strings
    @param hmm_database: The hmm database to use. Currently, all hmms are in the ALL database.
                         The code to develop new models is in build_pipeline in the git repo.
    @param ncpu: The number of cpu's to allow hmmer to use.
    """
    hmm_full_path = os.path.join(HMM_path, '{}.hmm'.format(hmm_database))
    if not os.path.exists(hmm_full_path):
        print("ANARCI HMM database not found:", hmm_full_path, file=sys.stderr)
        print("This indicates an error during installation. Please reinstall ANARCI.")
        sys.exit(1)

    with pyhmmer.plan7.HMMFile(hmm_full_path) as hmm_file:
        hmms = list(hmm_file)

    results = []
    hsps_per_sequence = _find_pyhmmer_hsps(sequence_list, hmms, ncpu=ncpu, bit_score_threshold=bit_score_threshold)
    assert len(sequence_list) == len(hsps_per_sequence), \
        'Unexpected mismatch in returned sequence hits: {} != {}'.format(len(sequence_list), len(hsps_per_sequence))

    for (sequence_id, sequence), hsps in zip(sequence_list, hsps_per_sequence):
        results.append(_parse_hmmer_query(
            hsps=hsps,
            seq_len=len(sequence),
            hmmer_species=hmmer_species
        ))
    return results


def _find_pyhmmer_hsps(sequence_list, hmms, ncpu=None, bit_score_threshold=80):
    alphabet = hmms[0].alphabet
    Z = len(hmms) # database size, used to calculate correct e-value

    sequences = [pyhmmer.easel.TextSequence(sequence=seq, accession=str(i).encode()).digitize(alphabet)
                 for i, (id, seq) in enumerate(sequence_list)]

    # We are using hmmsearch because hmmscan is not available in pyhmmer
    # It works the same except it returns a list of hits for each HMM model in the database
    # We are using the Z database size parameter to make sure that the evalue corresponds to hmmscan
    hits_per_hmm = pyhmmer.hmmsearch(
        hmms,
        sequences,
        cpus=0 if ncpu is None else ncpu,
        Z=Z
    )

    hsps_batch = [[] for i in range(len(sequence_list))]

    # Distribute hits by sequence and filter out hits with low bitscore
    for hmm_hits in hits_per_hmm:
        for hit in hmm_hits:
            for hsp in hit.domains:
                if hsp.score >= bit_score_threshold:
                    hsps_batch[int(hit.accession.decode())].append(hsp)

    return hsps_batch


def get_hmm_length( species, ctype ):
    '''
    Get the length of an hmm given a species and chain type. 
    This tells us how many non-insertion positions there could possibly be in a domain (127 or 128 positions under imgt)
    '''
    try:
        return len(list(all_germlines['J'][ctype][species].values())[0].rstrip('-'))
    except KeyError:
        return 128


def number_sequence_from_alignment(state_vector, sequence, scheme="imgt", chain_type=None):
    """
    Given you have an alignment. Give back the numbering
    
    @param state_vector: List of states from the hmm. Effectively these are imgt columns but CDR3 has not been redone. 
    @param sequence: The original sequence string or list.
    @param scheme: The numbering scheme to apply
    @param chain_type: The type of chain to apply numbering for. Some schemes do not require this (IMGT). Others (e.g. Chothia/Wolfguy) do.
    
    @return: A list of numbering identifier / amino acids tuples over the domain that has been numbered. The indices of the start (inclusive) and end point (exclusive) in the sequence for the numbering 
    """
    scheme=scheme.lower()
    if scheme == "imgt":
        return number_imgt(state_vector, sequence)
    elif scheme == "chothia":
        if chain_type == "H":
            return number_chothia_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_chothia_light(state_vector, sequence)
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))
    elif scheme == "kabat":
        if chain_type == "H":
            return number_kabat_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_kabat_light(state_vector, sequence)
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))
    elif scheme == "martin":
        if chain_type == "H":
            return number_martin_heavy(state_vector, sequence)
        elif chain_type in "KL":
            return number_martin_light(state_vector, sequence)
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))
    elif scheme == "aho":
        return number_aho(state_vector, sequence, chain_type) # requires the chain type to heuristically put the CDR1 gap in position.
    elif scheme == "wolfguy":
        if chain_type == "H":
            return number_wolfguy_heavy( state_vector, sequence )     
        elif chain_type in "KL":
            return number_wolfguy_light( state_vector, sequence )  
        else:
            raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))      
    else:
        raise AssertionError("Unimplemented numbering scheme %s for chain %s"%( scheme, chain_type))

def number_sequences_from_alignment(sequences, alignments, scheme="imgt", allow=set(["H","K","L","A","B","G","D"]), 
                                    assign_germline=False, allowed_species=None):
    '''
    Given a list of sequences and a corresponding list of alignments from run_hmmer apply a numbering scheme.
    '''

    # Iteration over the sequence alignments performing the desired numbering 
    numbered = []
    alignment_details = []
    hit_tables = []
    for i in range(len(sequences)):

        # Unpack
        hit_table, state_vectors, detailss = alignments[i] # We may have multiple domains per sequence (e.g. single chain fvs). 

        # Iterate over all the domains in the sequence that have been recognised (typcially only 1 with the current hmms available)
        hit_numbered, hit_details = [], []
        for di in range( len( state_vectors ) ):
            state_vector = state_vectors[di]
            details      = detailss[di]
            details["scheme"]=scheme
            details["query_name"]=sequences[i][0]            

            # Only number things that are allowed. We still keep the alignment details and hit_table
            if state_vector and details["chain_type"] in allow: 
                try:
                    # Do the numbering and validate (for development purposes)
                    hit_numbered.append( validate_numbering(number_sequence_from_alignment(state_vector, sequences[i][1], 
                                                            scheme=scheme, chain_type=details["chain_type"]), sequences[i] ) )
                    if assign_germline:
                        details["germlines"] = run_germline_assignment( state_vector, sequences[i][1], 
                                                                        details["chain_type"], allowed_species=allowed_species)
                    hit_details.append( details )
                except AssertionError as e: # Handle errors. Those I have implemented should be assertion.
                    print(str(e), file=sys.stderr)
                    raise e # Validation went wrong. Error message will go to stderr. Want this to be fatal during development.
                except Exception as e:
                    traceback.print_exc()
                    print("Error: Something really went wrong that has not been handled", file=sys.stderr)
                    print(str(e), file=sys.stderr)
                    raise e
                
        if hit_numbered: 
            numbered.append( hit_numbered )
            alignment_details.append( hit_details )
        else: 
            numbered.append( None )
            alignment_details.append( None )
        hit_tables.append(hit_table)

    return numbered, alignment_details, hit_tables

def get_identity( state_sequence, germline_sequence ):
    """
    Get the partially matched sequence identity between two aligned sequences. 
    Partial in the sense that gaps can be in the state_sequence.
    """
    # Ensure that the sequences are the expected length
    assert len( state_sequence) == len(germline_sequence ) == 128
    n, m = 0, 0
    for i in range( 128 ):
        if germline_sequence[i] == "-":continue
        if state_sequence[i].upper() == germline_sequence[i]: m+=1
        n+=1

    if not n:
        return 0    
    return float(m)/n
    

def run_germline_assignment(state_vector, sequence, chain_type, allowed_species=None ):
    """
    Find the closest sequence identity match.
    """
    genes={'v_gene': [None,None],
           'j_gene': [None,None],
         }


    # Extract the positions that correspond to match (germline) states. 
    state_dict = dict( ((i, 'm'),None) for i in range(1,129))
    state_dict.update(dict(state_vector))
    state_sequence = "".join([ sequence[state_dict[(i, 'm')]] if state_dict[(i,'m')] is not None else "-" for i in range(1,129) ])

    # Iterate over the v-germline sequences of the chain type of interest.
    # The maximum sequence identity is used to assign the germline 
    if chain_type in all_germlines["V"]:
        if allowed_species is not None:
            if not all( [ sp in all_germlines['V'][chain_type] for sp in allowed_species ] ): # Made non-fatal
                return {}
        else:
            allowed_species = all_species
        seq_ids = {}
        for species in allowed_species:
            if species not in all_germlines["V"][ chain_type ]: continue # Previously bug.
            for gene, germline_sequence in all_germlines["V"][ chain_type ][ species ].items():
                seq_ids[ (species, gene) ] = get_identity( state_sequence , germline_sequence )
        genes['v_gene' ][0] = max( seq_ids, key=lambda x: seq_ids[x] )
        genes['v_gene' ][1] = seq_ids[ genes['v_gene' ][0] ]
        
        # Use the assigned species for the v-gene for the j-gene. 
        # This assumption may affect exotically engineered abs but in general is fair.
        species = genes['v_gene' ][0][0]       
        if chain_type in all_germlines["J"]:
            if species in all_germlines["J"][chain_type]:
                seq_ids = {}
                for gene, germline_sequence in all_germlines["J"][ chain_type ][ species ].items():
                    seq_ids[ (species, gene) ] = get_identity( state_sequence , germline_sequence )
                genes['j_gene' ][0] = max( seq_ids, key=lambda x: seq_ids[x] )
                genes['j_gene' ][1] = seq_ids[ genes['j_gene' ][0] ]
     
    return genes

def check_for_j( sequences, alignments, scheme ):
    '''
    As the length of CDR3 gets long (over 30ish) an alignment that does not include the J region becomes more favourable.
    This leads to really long CDR3s not being numberable. 

    To overcome this problem, when no J region is detected we try without the v region.
    '''
    for i in range( len( sequences ) ):
        # Check the alignment for J region
        if len(alignments[i][1]) ==1: # Only do for single domain chains. 

            # Check whether a J region has been identified. If not check whether there is still a considerable amount of sequence   
            # remaining. 
            ali = alignments[i][1][0]

            # Find the last match position. 
            last_state  = ali[-1][0][0]
            last_si     = ali[-1][1]
            if last_state < 120: # No or very little J region
                if last_si + 30 < len( sequences[i][1] ): # Considerable amount of sequence left...suspicious of a long CDR3
                    # Find the position of the conserved cysteine (imgt 104). 
                    cys_si = dict( ali ).get( (104,'m'), None )
                    if cys_si is not None: # 104 found.

                        # Find the corresponding index in the alignment.
                        cys_ai = ali.index( ((104, 'm'), cys_si) )
                
                        # Try to identify a J region in the remaining sequence after the 104. A low bit score threshold is used.
                        _, re_states, re_details  = run_pyhmmer( [(sequences[i][0], sequences[i][1][cys_si+1:])],
                                                               bit_score_threshold=10 )[0] 

                        # Check if a J region was detected in the remaining sequence.
                        if re_states and re_states[0][-1][0][0] >= 126 and re_states[0][0][0][0] <= 117: 

                            # Sandwich the presumed CDR3 region between the V and J regions.

                            vRegion   = ali[:cys_ai+1]
                            jRegion   = [ (state, index+cys_si+1) for state, index in re_states[0] if state[0] >= 117 ]
                            cdrRegion = []
                            next = 105
                            for si in range( cys_si+1, jRegion[0][1] ):
                                if next >= 116:
                                    cdrRegion.append( ( (116, 'i'), si ) )
                                else:
                                    cdrRegion.append( ( (next, 'm'), si ) )
                                    next +=1 

                            # Update the alignment entry.
                            alignments[i][1][0] = vRegion + cdrRegion + jRegion 
                            alignments[i][2][0]['query_end'] = jRegion[-1][1] + 1 



##################################
# High level numbering functions #
##################################

# Main function for ANARCI 
# Name conflict with function, module and package is kept for legacy unless issues are reported in future. 
def anarci(sequences, scheme="imgt", database="ALL", output=False, outfile=None, csv=False, allow=set(["H","K","L","A","B","G","D"]), 
           ncpu=None, assign_germline=False, allowed_species=['human','mouse'], bit_score_threshold=80):
    """
    The main function for anarci. Identify antibody and TCR domains, number them and annotate their germline and species. 

    It is advised to use one of the wrapper functions:
        o run_anarci   - fasta file or sequence list in. Automated multiprocessing for large jobs. Sequences, numbering, details 
                         and hit tables out. 
        o number       - single sequence in, numbering out

    
    @param sequences: A list or tuple of (Id, Sequence) pairs
                              e.g. [ ("seq1","EVQLQQSGAEVVRSG ..."),
                                     ("seq2","DIVMTQSQKFMSTSV ...") ]
    @param scheme:    The numbering scheme that should be applied. Choose from imgt, chothia, kabat or martin
    @param output:    Boolean flag to say whether the result should be output.
    @param outfile:   The name of the file to output to. If output is True and outfile is None then output is printed
                      to stdout.
    @param csv:       Boolean flag to say whether the csv output alignment format or the vertical anarci format should be used.
    @param allow:     A set containing the chain types that should be recognised. If chothia, kabat or martin is used
                      as the scheme, anarci will ignore tcr chains. Choose a subset of ["H","K","L","A","B","G","D"]
    @param assign_germline: Using highest sequence identity assign the germline to the chain. Can be more accurate at identifying
                      species than the best HMM hit alone. (Bool)
    @param allowed_species: If assign_germline is true, limit the species that can be assigned to a limited set. Useful when the 
                      animal species is known or when performing closest germline experiments. Choose a subset of ['human',
                      'mouse','rat','rabbit','rhesus','pig','alpaca'].


    @param bit_score_threshold: The threshold score from HMMER at which an alignment should be numbered. Lowering the threshold 
                      means domain recognition is more permissive and can be useful for numbering heavily engineered molecules. 
                      However, too low and false positive recognition of other ig-like molecules will occur.
    @param ncpu:      The number of cpu's that hmmer should be allowed to use. If not specified then the hmmscan
                      default is used. N.B. hmmscan must be compiled with multithreading enabled for this option to have effect. 
                      Please consider using the run_anarci function for native multiprocessing with anarci.
    @param database:  The HMMER database that should be used. Normally not changed unless a custom db is created.


    @return: Three lists. Numbered, Alignment_details and Hit_tables.
             Each list is in the same order as the input sequences list.
             A description of each entry in the three lists is as followed.
               o Numbered: will be None if no domain was found for that sequence or a list of domains with their 
                           numbering, start and finish indices.
               o Alignment_details: will be None if no domain was found for that sequence or a dictionary for each
                           domain identified containing the details of the alignment (chain type, e-value, species etc).
               o Hit_tables: None if no domain was found for that sequence or a nested list for each domain containing
                           the hit table from hmmscan.
    
    """
    
    # Validate the input scheme
    try:
        scheme = scheme_short_to_long[scheme.lower()]
    except KeyError:
        raise AssertionError("Unrecognised or unimplemented scheme: %s"%scheme)        

    # Check we have arguments for output before doing work.
    if csv:
        assert outfile, 'If csv output is True then an outfile must be specified'
        _path, _ = os.path.split(outfile)
        assert (not _path) or os.path.exists(_path), 'Output directory %s does not exist'%_path


    # Perform the alignments of the sequences to the hmm database
    alignments = run_pyhmmer(sequences,hmm_database=database,ncpu=ncpu,bit_score_threshold=bit_score_threshold,hmmer_species=allowed_species )
     
    # Check the numbering for likely very long CDR3s that will have been missed by the first pass.
    # Modify alignments in-place
    check_for_j( sequences, alignments, scheme )

    # Apply the desired numbering scheme to all sequences
    numbered, alignment_details, hit_tables = number_sequences_from_alignment(sequences, alignments, scheme=scheme, allow=allow, 
                                                                              assign_germline=assign_germline, 
                                                                              allowed_species=allowed_species)

    # Output if necessary
    if output: 
        if csv:
            csv_output(sequences, numbered, alignment_details, outfile)
        else:
            outto, close=sys.stdout, False
            if outfile:
                outto, close = open(outfile,'w'), True
            anarci_output(numbered, sequences, alignment_details, outto)
            if close:
                outto.close()


    return numbered, alignment_details, hit_tables

# Wrapper to run anarci using multiple processes and automate fasta file reading.
def run_anarci( seq, ncpu=1, **kwargs):
    '''
    Run the anarci numbering protocol for single or multiple sequences.
    
    @param sequences: A list or tuple of (Id, Sequence) pairs
                              e.g. [ ("seq1","EVQLQQSGAEVVRSG ..."),
                                     ("seq2","DIVMTQSQKFMSTSV ...") ]
    @param scheme:    The numbering scheme that should be applied. Choose from imgt, chothia, kabat or martin
    @param output:    Boolean flag to say whether the result should be output.
    @param outfile:   The name of the file to output to. If output is True and outfile is None then output is printed
                      to stdout.
    @param allow:     A set containing the chain types that should be recognised. If chothia, kabat or martin is used
                      as the scheme, anarci will ignore tcr chains. Choose a subset of ["H","K","L","A","B","G","D"]
    @param assign_germline: Using highest sequence identity assign the germline to the chain. Can be more accurate at identifying
                      species than the best HMM hit alone. (Bool)
    @param allowed_species: If assign_germline is true, limit the species that can be assigned to a limited set. Useful when the 
                      animal species is known or when performing closest germline experiments. Choose a subset of ['human',
                      'mouse','rat','rabbit','rhesus','pig','alpaca'].

    @param bit_score_threshold: The threshold score from HMMER at which an alignment should be numbered. Lowering the threshold 
                      means domain recognition is more permissive and can be useful for numbering heavily engineered molecules. 
                      However, too low and false positive recognition of other ig-like molecules will occur.
    @param hmmerpath: The path to hmmscan. If left unspecified then the PATH will be searched. 
    @param ncpu:      The number of cpu's that hmmer should be allowed to use. If not specified then the hmmscan 
                      default is used. N.B. hmmscan must be compiled with multithreading enabled for this option to have effect. 
                      Please consider using the run_anarci function for native multiprocessing with anarci.
    @param database:  The HMMER database that should be used. Normally not changed unless a custom db is created.

    @return: Four lists. Sequences, Numbered, Alignment_details and Hit_tables.
             Each list is in the same order. 
             A description of each entry in the four lists is as followed.
               o Sequences: The list of sequences formatted as [(Id,sequence), ...]. 
               o Numbered: will be None if no domain was found for that sequence or a list of domains with their 
                           numbering, start and finish indices.
               o Alignment_details: will be None if no domain was found for that sequence or a dictionary for each
                           domain identified containing the details of the alignment (chain type, e-value, species etc).
               o Hit_tables: None if no domain was found for that sequence or a nested list for each domain containing
                           the hit table from hmmscan.

    '''
    # Parse the input sequence or fasta file.
    if isinstance(seq, list) or isinstance(seq,tuple): # A list (or tuple) of (name,sequence) sequences
        assert all( len(_) == 2 for _ in seq ), "If list or tuple supplied as input format must be [ ('ID1','seq1'), ('ID2', 'seq2'), ... ]"
        sequences = seq
    elif os.path.isfile( seq ): # Fasta file.
        # Read the sequences. All are read into memory currently...
        sequences = read_fasta( seq ) 
        ncpu = int(max(1, ncpu )) 
    elif isinstance(seq, str): # Single sequence
        validate_sequence( seq )
        ncpu=1
        sequences = [ ["Input sequence", seq ]]

    # Handle the arguments to anarci.
    output  = kwargs.get('output', False )
    outfile = kwargs.get('outfile', False )
    csv = kwargs.get( 'csv', False )
    if csv: # Check output arguments before doing work.
        assert outfile, 'If csv output is True then an outfile must be specified'
        _path, _ = os.path.split(outfile)
        assert (not _path) or os.path.exists(_path), 'Output directory %s does not exist'%_path
        
    kwargs['ncpu'] = 1 # Set hmmscan ncpu to 1. HMMER has to be compiled appropriately for this to have an effect. 
    kwargs['output'] = False # Overide and write the compiled results here. 

    anarci_partial = partial( anarci, **kwargs )        
    chunksize = math.ceil( float( len(sequences) )/ncpu )

    # Run the anarci function using a pool of workers. Using the map_async to get over the KeyboardInterrupt bug in python2.7
    if ncpu > 1:
        pool = Pool( ncpu )
        results = pool.map_async( anarci_partial, grouper( chunksize, sequences ) ).get()
        pool.close()
    else:
        results = list(map( anarci_partial, grouper( chunksize, sequences ) ))

    # Reformat the results to flat lists.
    numbered = sum( (_[0] for _ in results), [] )
    alignment_details = sum( (_[1] for _ in results ), [] )
    hit_tables = sum( (_[2] for _ in results), [] )

    # Output if necessary
    if output: 
        if csv:
            csv_output(sequences, numbered, alignment_details, outfile)             
        else:
            outto, close=sys.stdout, False
            if outfile:
                outto, close = open(outfile,'w'), True
            anarci_output(numbered, sequences, alignment_details, outto)
            if close:
                outto.close()

    # Return the results
    return sequences, numbered, alignment_details, hit_tables
                


# Wrapper function for simple sequence in numbering and chain type out behaviour. 
def number(sequence, scheme="imgt", database="ALL", allow=set(["H","K","L","A","B","G","D"]), allowed_species=['human','mouse']):
    """
    Given a sequence string, use anarci to number it using the scheme of choice.
    Only the first domain will be recognised and numbered

    For multiple sequences it is advised to use run_anarci instead of iterative use of this function.

    @param sequence: An amino acid sequence string
    @param scheme: The numbering scheme that should be applied. Choose from imgt, chothia, kabat or martin
    @param database: The HMMER database that should be used. Normally not changed unless a custom db is created.
    @param allow: A set containing the chain types that should be recognised. If chothia, kabat or martin is used
                  as the scheme, anarci will ignore tcr chains.

    @return: If the sequence can be numbered, a list containing the numbering and sequence; and the chain type. 
             Otherwise both are False.

    """
    
    try:
        validate_sequence( sequence )  
        scheme = scheme_short_to_long[scheme.lower()]
    except KeyError:
        raise AssertionError("Unrecognised to unimplemented scheme: %s"%scheme)
    
    if len(sequence) < 70: # Length check. ANARCI can number fragments of chains well. Encourage full domain numbering. 
        return False, False
   
    try:
        if not allowed_species:
            numbered, alignment_details, _ = anarci( [("sequence_0", sequence)], scheme=scheme, database=database, output=False, allow=allow )
        else:
            numbered, alignment_details, _ = anarci( [("sequence_0", sequence)], scheme=scheme, database=database, output=False, allow=allow, allowed_species = allowed_species )
    except AssertionError: # Catch where the user has tried to number a TCR with an antibody scheme
        return False, False
    

    # We return the numbering list and the chain type where kappa and lambda chains are both "L" for light
    if numbered[0]:
        return numbered[0][0][0], chain_type_to_class[alignment_details[0][0]["chain_type"]]
    else:
        return False, False

if __name__ == "__main__":
    # Test and example useage of the anarci function. 
    sequences = [ ("12e8:H","EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAP"),
                  ("12e8:L","DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASV"),
                  ("scfv:A","DIQMTQSPSSLSASVGDRVTITCRTSGNIHNYLTWYQQKPGKAPQLLIYNAKTLADGVPSRFSGSGSGTQFTLTISSLQPEDFANYYCQHFWSLPFTFGQGTKVEIKRTGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKRLEWVAYISSGGGSTYFPDTVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARQNKKLTWFDYWGQGTLVTVSSHHHHHH"),
                  ("lysozyme:A","KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL")]

    results = anarci(sequences, scheme="imgt", output=True)
    numbering, alignment_details, hit_tables = results

    expect_one_VH_domain_numbering, expect_one_VL_domain_numbering, expect_VH_then_VL_numbering, expect_None = numbering
    assert  len(expect_one_VH_domain_numbering) == 1
    assert  len(expect_one_VL_domain_numbering) == 1
    assert  len(expect_VH_then_VL_numbering)    == 2
    assert  expect_None                         == None




