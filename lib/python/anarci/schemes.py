#    ANARCI - Antibody Numbering and Antigen Receptor ClassIfication
#    Copyright (C) 2016 Oxford Protein Informatics Group (OPIG)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the BSD 3-Clause License.
#
#    You should have received a copy of the BSD 3-Clause Licence
#    along with this program.  If not, see <https://opensource.org/license/bsd-3-clause/>.

'''
Module containing functions to convert hmm alignment to a numbering scheme. 

Currently implemented

For IG's
IMGT
Chothia
Kabat
Martin (Extended Chothia)
Aho 
Wolfguy

For TR's
IMGT
(Aho)

---------------------------------------------------------------------------------------------------------------------
Functions are written to a template:

There are 128 match states in the HMMs (these are the IMGT states). The alignment to these states must be converted to
correspond to the scheme of choice. 

We define:
  - a state string consisting of 'X' and 'I' where:
    X  means that for the state there is an equivalent position in the numbering scheme.
    I  means that for the state there is not an equivalent position in the numbering scheme. It should therefore be 
       considered as an insertion in the scheme.
       
  - a region string consisting of characters (integers in the currently implemented schemes). Each character 
corresponds to a contiguous region. Therefore each state can be assigned a region according to the scheme. 

  - a mapping between region characters and region indices as a dictionary. e.g. the first region character maps
to 0, second to 1 ...

  - a dictionary containing the difference between state number (imgt) and scheme number at the *beginning* of 
each region using the region indices as keys and the difference as values.

  - the number of regions defined

  - a list for which delete states should not be included in the numbering (typically those for the cdrs). This
will allow the length of the region to be the number of residues found instead of the number of possible states plus
insertions.

 
This all goes into the _number_regions function along with the sequence and the state_vector (the alignment from the
HMM).

_number regions will then divide the aligned part of the sequence into as many regions as defined above. Within each 
region it will give a numbering according to the input parameters. A list of lists will be returned containing the 
numbered sequence for each region.

Some of the regions will not be numbered correctly according to the scheme. For example the insertions for the CDRs
will not necessarily be on the correct residue. For each different scheme these regions are then modified (see code
for implementation) 

Finally the full numbered sequence is compiled and returned to the calling function.
---------------------------------------------------------------------------------------------------------------------

Other schemes can be implemented following the template above. 


'''

# Alphabet used for insertion (last (-1th) is a blank space for no insertion)
alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ", "KK", "LL", "MM", "NN", "OO", "PP", "QQ", "RR", "SS", "TT", "UU", "VV", "WW", "XX", "YY", "ZZ", " "]

# Blosum62 matrix. Used in some annotation methods to recognise pre-defined motifs
blosum62 = {('B', 'N'): 3, ('W', 'L'): -2, ('G', 'G'): 6, ('X', 'S'): 0, ('X', 'D'): -1, ('K', 'G'): -2, ('S', 'E'): 0, ('X', 'M'): -1, ('Y', 'E'): -2, ('W', 'R'): -3, ('I', 'R'): -3, ('X', 'Z'): -1, ('H', 'E'): 0, ('V', 'M'): 1, ('N', 'R'): 0, ('I', 'D'): -3, ('F', 'D'): -3, ('W', 'C'): -2, ('N', 'A'): -2, ('W', 'Q'): -2, ('L', 'Q'): -2, ('S', 'N'): 1, ('Z', 'K'): 1, ('V', 'N'): -3, ('Q', 'N'): 0, ('M', 'K'): -1, ('V', 'H'): -3, ('G', 'E'): -2, ('S', 'L'): -2, ('P', 'R'): -2, ('D', 'A'): -2, ('S', 'C'): -1, ('E', 'D'): 2, ('Y', 'G'): -3, ('W', 'P'): -4, ('X', 'X'): -1, ('Z', 'L'): -3, ('Q', 'A'): -1, ('V', 'Y'): -1, ('W', 'A'): -3, ('G', 'D'): -1, ('X', 'P'): -2, ('K', 'D'): -1, ('T', 'N'): 0, ('Y', 'F'): 3, ('W', 'W'): 11, ('Z', 'M'): -1, ('L', 'D'): -4, ('M', 'R'): -1, ('Y', 'K'): -2, ('F', 'E'): -3, ('M', 'E'): -2, ('S', 'S'): 4, ('X', 'C'): -2, ('Y', 'L'): -1, ('H', 'R'): 0, ('P', 'P'): 7, ('K', 'C'): -3, ('S', 'A'): 1, ('P', 'I'): -3, ('Q', 'Q'): 5, ('L', 'I'): 2, ('P', 'F'): -4, ('B', 'A'): -2, ('Z', 'N'): 0, ('M', 'Q'): 0, ('V', 'I'): 3, ('Q', 'C'): -3, ('I', 'H'): -3, ('Z', 'D'): 1, ('Z', 'P'): -1, ('Y', 'W'): 2, ('T', 'G'): -2, ('B', 'P'): -2, ('P', 'A'): -1, ('C', 'D'): -3, ('Y', 'H'): 2, ('X', 'V'): -1, ('B', 'B'): 4, ('Z', 'F'): -3, ('M', 'L'): 2, ('F', 'G'): -3, ('S', 'M'): -1, ('M', 'G'): -3, ('Z', 'Q'): 3, ('S', 'Q'): 0, ('X', 'A'): 0, ('V', 'T'): 0, ('W', 'F'): 1, ('S', 'H'): -1, ('X', 'N'): -1, ('B', 'Q'): 0, ('K', 'A'): -1, ('I', 'Q'): -3, ('X', 'W'): -2, ('N', 'N'): 6, ('W', 'T'): -2, ('P', 'D'): -1, ('B', 'C'): -3, ('I', 'C'): -1, ('V', 'K'): -2, ('X', 'Y'): -1, ('K', 'R'): 2, ('Z', 'R'): 0, ('W', 'E'): -3, ('T', 'E'): -1, ('B', 'R'): -1, ('L', 'R'): -2, ('Q', 'R'): 1, ('X', 'F'): -1, ('T', 'S'): 1, ('B', 'D'): 4, ('Z', 'A'): -1, ('M', 'N'): -2, ('V', 'D'): -3, ('F', 'A'): -2, ('X', 'E'): -1, ('F', 'H'): -1, ('M', 'A'): -1, ('K', 'Q'): 1, ('Z', 'S'): 0, ('X', 'G'): -1, ('V', 'V'): 4, ('W', 'D'): -4, ('X', 'H'): -1, ('S', 'F'): -2, ('X', 'L'): -1, ('B', 'S'): 0, ('S', 'G'): 0, ('P', 'M'): -2, ('Y', 'M'): -1, ('H', 'D'): -1, ('B', 'E'): 1, ('Z', 'B'): 1, ('I', 'E'): -3, ('V', 'E'): -2, ('X', 'T'): 0, ('X', 'R'): -1, ('R', 'R'): 5, ('Z', 'T'): -1, ('Y', 'D'): -3, ('V', 'W'): -3, ('F', 'L'): 0, ('T', 'C'): -1, ('X', 'Q'): -1, ('B', 'T'): -1, ('K', 'N'): 0, ('T', 'H'): -2, ('Y', 'I'): -1, ('F', 'Q'): -3, ('T', 'I'): -1, ('T', 'Q'): -1, ('P', 'L'): -3, ('R', 'A'): -1, ('B', 'F'): -3, ('Z', 'C'): -3, ('M', 'H'): -2, ('V', 'F'): -1, ('F', 'C'): -2, ('L', 'L'): 4, ('M', 'C'): -1, ('C', 'R'): -3, ('D', 'D'): 6, ('E', 'R'): 0, ('V', 'P'): -2, ('S', 'D'): 0, ('E', 'E'): 5, ('W', 'G'): -2, ('P', 'C'): -3, ('F', 'R'): -3, ('B', 'G'): -1, ('C', 'C'): 9, ('I', 'G'): -4, ('V', 'G'): -3, ('W', 'K'): -3, ('G', 'N'): 0, ('I', 'N'): -3, ('Z', 'V'): -2, ('A', 'A'): 4, ('V', 'Q'): -2, ('F', 'K'): -3, ('T', 'A'): 0, ('B', 'V'): -3, ('K', 'L'): -2, ('L', 'N'): -3, ('Y', 'N'): -2, ('F', 'F'): 6, ('L', 'G'): -4, ('B', 'H'): 0, ('Z', 'E'): 4, ('Q', 'D'): 0, ('X', 'B'): -1, ('Z', 'W'): -3, ('S', 'K'): 0, ('X', 'K'): -1, ('V', 'R'): -3, ('K', 'E'): 1, ('I', 'A'): -1, ('P', 'H'): -2, ('B', 'W'): -4, ('K', 'K'): 5, ('H', 'C'): -3, ('E', 'N'): 0, ('Y', 'Q'): -1, ('H', 'H'): 8, ('B', 'I'): -3, ('C', 'A'): 0, ('I', 'I'): 4, ('V', 'A'): 0, ('W', 'I'): -3, ('T', 'F'): -2, ('V', 'S'): -2, ('T', 'T'): 5, ('F', 'M'): 0, ('L', 'E'): -3, ('M', 'M'): 5, ('Z', 'G'): -2, ('D', 'R'): -2, ('M', 'D'): -3, ('W', 'H'): -2, ('G', 'C'): -3, ('S', 'R'): -1, ('S', 'I'): -2, ('P', 'Q'): -1, ('Y', 'A'): -2, ('X', 'I'): -1, ('E', 'A'): -1, ('B', 'Y'): -3, ('K', 'I'): -3, ('H', 'A'): -2, ('P', 'G'): -2, ('F', 'N'): -3, ('H', 'N'): 1, ('B', 'K'): 0, ('V', 'C'): -1, ('T', 'L'): -1, ('P', 'K'): -1, ('W', 'S'): -3, ('T', 'D'): -1, ('T', 'M'): -1, ('P', 'N'): -2, ('K', 'H'): -1, ('T', 'R'): -1, ('Y', 'R'): -2, ('L', 'C'): -1, ('B', 'L'): -4, ('Z', 'Y'): -2, ('W', 'N'): -4, ('G', 'A'): 0, ('S', 'P'): -1, ('E', 'Q'): 2, ('C', 'N'): -3, ('H', 'Q'): 0, ('D', 'N'): 1, ('Y', 'C'): -2, ('L', 'H'): -3, ('E', 'C'): -4, ('Z', 'H'): 0, ('H', 'G'): -2, ('P', 'E'): -1, ('Y', 'S'): -2, ('G', 'R'): -2, ('B', 'M'): -3, ('Z', 'Z'): 4, ('W', 'M'): -1, ('Y', 'T'): -2, ('Y', 'P'): -3, ('Y', 'Y'): 7, ('T', 'K'): -1, ('Z', 'I'): -3, ('T', 'P'): -1, ('V', 'L'): 1, ('F', 'I'): 0, ('G', 'Q'): -2, ('L', 'A'): -1, ('M', 'I'): 1}


def smooth_insertions(state_vector):
    '''
    The function aims to correct to the expected imgt alignment. Renumbering functions then translate from the imgt scheme to the
    appropriate scheme.        

    Handle insertions made by HMMER that we suspect may be in the wrong position.
    Edge cases include:
        - Insertions at the C terminal of fw1, fw3 and fw3 regions. Can occur when 'conserved' residues have been mutated and the 
          same amino acid appears in the the following CDR (e.g. mutate cysteine at 104 but the CDR3 has one or more cysteines)
        - Same as above possible (but not observed in structure seqs) for N terminal of fw2, fw3 and fw4... TODO
        - Heavily mutated N terminal regions that are partially recognised (e.g. 3gk8 chain H). Insertions should not be allowed 
          before N terminal deletions have been used. Preserve deletion locations that are not N terminal (e.g. 10 in IMGT H) if 
          the gap has been placed by the alignment.

    '''
    # Small overhead doing these corrections but worth it for reducing edge cases.
    
    # Enforce insertion patterns as below. The CDRs are renumbered in each case so that insertions are placed accoring to the scheme    
#  '11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777'
#  '                        mmmi                         mmmi                                             mmmi                      '
#  '                        mmmi        immm             mmmi      immm                                   mmmi         immm         '

    # Enforce any insertions at the end and beginning of framework regions to be moved into the CDR region for renumbering. 
    enforced_patterns = [ [(25,'m'),(26,'m'),( 27,'m'),( 28,'i')],
                          [(38,'i'),(38,'m'),(39,'m'),(40,'m')],
                          [(54,'m'),(55,'m'),(56,'m'),(57,'i')],
                          [(65,'i'),(65,'m'),(66,'m'),(67,'m')],
                          [(103,'m'),(104,'m'),(105,'m'),(106,'i')],
                          [(117,'i'),(117,'m'),(118,'m'),(119,'m')] ]

    # Insertions in FW1 are only allowed if there are a fewer number of n-terminal deletions made. 

    state_buffer = [] 
    sv = []
    for (state_id, state_type ), si in state_vector:
        if state_id < 23: # Everything before the cysteine at 23.
            state_buffer.append( ((state_id, state_type ), si) )
            reg = -1   
        elif 25 <= state_id < 28: # Add to the buffer 
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 0
        elif 37 < state_id <= 40: # Add to the buffer 
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 1
        elif 54 <= state_id < 57: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 2
        elif 64 < state_id <= 67: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 3
        elif 103 <= state_id < 106: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 4
        elif 116 < state_id <= 119: # Add to the buffer
            state_buffer.append( ((state_id, state_type ), si) )
            reg = 5
        elif len(state_buffer) != 0: # Add the buffer and reset

            # Find the number of insertions in the buffer
            nins = sum( 1 for s in state_buffer if s[0][1] == 'i' ) 

            # If there are insertions, adjust the alignment
            if nins > 0: # We have insertions

                if reg == -1: # FW1, only adjust if there are the same or more N terminal deletions than insertions
                    nt_dels = state_buffer[0][0][0] - 1 # Missing states
                    for (_id, _type ), _si in state_buffer: # Explicit deletion states.
                        if _type == 'd' or _si == None:
                            nt_dels +=1 
                        else: # First residue found
                            break
                    if nt_dels >= nins: # More n terminal deletions than insertions found. Likely misalignment.
                        
                        # Preserve the deleted states structure by using the same match annotations
                        new_states = [ s for s, _ in state_buffer if s[1] == 'm'] 
                        _first = new_states[0][0]

                        # Remove the deletions so that only residue positions are included
                        state_buffer = [ s for s in state_buffer if s[0][1] != 'd' ]

                        # Extend N terminal states backwards from the first match states
                        _add = len( state_buffer ) - len( new_states ) 
                        assert _add >= 0, 'Implementation logic error' # Should be adding a positive number of positions
                        new_states = [ (_,'m') for _ in range( _first - _add, _first ) ] + new_states
                        assert len(new_states)==len(state_buffer), 'Implementation logic error' # Should have the same length

                        # Assign them preserving the order of the sequence. 
                        for i in range( len(state_buffer ) ):
                            sv.append( ( new_states[i], state_buffer[i][1]) )
                    else:
                        sv += state_buffer # The insertions may be incorrect but unknown what to do. Let the alignment place.
                else:
                    # Remove any deletions in the buffer. Unlikely to happen but do anyway
                    state_buffer = [ s for s in state_buffer if s[0][1] != 'd' ]
        
                    # Define the new states defined by the enforced pattern and the length of the buffer
                    if reg % 2: # nterm fw
                        new_states = [enforced_patterns[reg][0]]*max( 0, len(state_buffer)-3) + enforced_patterns[reg][ max( 4-len(state_buffer), 1):]
                    else: # cterm fw
                        new_states = enforced_patterns[reg][:3] + [enforced_patterns[reg][2]]*max( 0, len(state_buffer)-3)
                    # Assign them preserving the order of the sequence. 
                    for i in range( len(state_buffer ) ):
                        sv.append( ( new_states[i], state_buffer[i][1]) )
                                
            else: # Nothing to do - either all match or deletion states.
                sv += state_buffer

            # Add the current state
            sv.append( ((state_id, state_type ), si) )

            # Reset state buffer
            state_buffer = [] 
            
        else: # Simply append 
            sv.append( ((state_id, state_type ), si) )
    

    return sv


# General function to give annotations for regions that have direct mappings onto the hmm alignment (imgt states)
def _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions):
    """
    General function to number a sequence and divide it into different regions  
    
    @param sequence: The sequence string
    @param state_vector: The list of states from the aligned hmm
    @param state_string: A string of states for the scheme relative to IMGT (this is X for a direct equivalence, I if needs to be treated as insertion)
    @param region_string: A string of characters that indicate which hmm states are in each regions for this scheme (i.e. how should the sequence be divided up)
    @param region_index_dict: A dictionary converting the characters in region string to an index of the regions. 
    @param rels: The difference of the numbering integer at the *start* of each region
    @param n_regions: The number of regions
    @param exclude_deletions: A list of region indices for which deletion states should not be included. Typically the CDRs. 
                              These will be reannotated in the scheme function. Also allows the reset of insertions. 
    
    @return: A list of lists where each region has been numbered according to the scheme. Some regions will need renumbering. This should be taken care of after the function called.
    
    """

    state_vector = smooth_insertions( state_vector )
    
    _regions = [ [] for _ in range(n_regions) ]
    
    # Initialise the insertion index (-1 is a blank space) and the previous state.
    insertion = -1
    previous_state_id = 1
    previous_state_type = 'd'
    start_index, end_index  = None, None
    
    region = None

    # Iterate over the aligned state vector
    for (state_id, state_type ), si in state_vector:
       
        # Retrieve the region index
        if state_type != "i" or region is None: # BUG_FIX - JD 9/4/15 - do not allow a new region to start as an insertion.
            region = region_index_dict[region_string[state_id-1]] 

       
        # Check the state_types
        if state_type == "m": # It is a match
            
            # Check whether this position is in the scheme as an independent state
            if state_string[state_id-1]=="I": # No, it should be treated as an insertion
                if previous_state_type != 'd': # Unless there was a deletion beforehand in which case this should be a real pos.
                    insertion +=1 # Increment the insertion annotation index
                rels[region] -= 1 # Update the relative numbering from the imgt states
            else: # Yes 
                insertion = -1 # Reset the insertions 
            
            # Add the numbering annotation to the appropriate region list            
            _regions[region].append( ( (state_id + rels[region], alphabet[insertion] ), sequence[si]  ) )
            previous_state_id = state_id # Record the previous state ID
            if start_index is None:
                start_index = si
            end_index = si

            previous_state_type = state_type
            
        elif state_type == "i": # It is an insertion
            insertion +=1 # Increment the insertion annotation index
            
            # Add the numbering annotation to the appropriate region list
            _regions[region].append( ( (previous_state_id + rels[region], alphabet[insertion]), sequence[si]  ) )
            if start_index is None:
                start_index = si
            end_index = si

            previous_state_type = state_type

        else: # It is a deletion
            previous_state_type = state_type
            
            # Check whether this position is in the scheme as an independent state
            if state_string[state_id-1]=="I": # No, therefore irrelevant to the scheme.
                rels[region] -= 1 # Update the relative numbering from the imgt states
                continue 
            
            insertion = -1 # Reset the insertions
            previous_state_id = state_id # Record the previous state ID, should not be needed (no delete to insert state transition)


        # Reset the inssertion index if necessary and allowed. (Means the insertion code is meaningless and will be reannotated)
        if insertion >= 25 and region in exclude_deletions:
            insertion = 0 
        
        assert insertion < 25, "Too many insertions for numbering scheme to handle" # We ran out of letters.
            
    return _regions, start_index, end_index


# Functions to perform the numbering and the corrections for each of the implemented schemes.
# These have been written fairly verbosely so that the template of how to generate a function for a new scheme is more clear.
# They have two stages: Perform the mapping between imgt and the scheme; Renumber those regions that do not map nicely onto imgt (e.g. CDR insertions)
    

    
########
# IMGT # 
########
# - Renumbering of the CDR 1 and 2 regions in IMGT has now been implemented to ensure consistency with the gapping rules of the 
# scheme. Previously gaps were defined using the HMM alignment as the underlying model was already based on the IMGT scheme. This 
# worked well in original test cases but appears to give inaccurate annotations in a significant number of cases in NGS size 
# sequence sets. We therefore now explicitly renumber the CDR 1 and 2 as with all the other schemes. 

def number_imgt(state_vector, sequence):
    """    
    Apply the IMGT numbering scheme for heavy or light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in IMGT scheme, I is an insertion. (All X's for IMGT)
    XXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXX XXXXXXXXXXXXXXXXX XXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXX
    11111111111111111111111111 222222222222 33333333333333333 4444444444 555555555555555555555555555555555555555 6666666666666 77777777777

    Regions - (N.B These do not match up with any particular definition of CDR)
    1. All positions before CDR1
    2. CDR1 positions
    3. Positions between CDR1/2
    4. CDR2 positions
    5. Positions between CDR2/3
    6. CDR positions 105 (inc) to 118 (exc)
    7. Positions after CDR3    
    
    """
    
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777'

    region_index_dict = {
                         "1":0,
                         "2":1,
                         "3":2,
                         "4":3,
                         "5":4,
                         "6":5,
                         "7":6
                         }
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                          1:0,
                          2:0,
                          3:0,
                          4:0,
                          5:0,
                          6:0,
                          7:0
                          }
    
    n_regions = 7

    exclude_deletions = [1,3,5]    
    
    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############

    _numbering = [ _regions[0], # Fw1
                   [],          # CDR1
                   _regions[2], # Fw2
                   [],          # CDR2
                   _regions[4], # Fw3
                   [],          # CDR3
                   _regions[6], # Fw4

                 ]

    # The alignment from HMMER should be correct for CDRs 1 and 2. Testing has shown not always the case and 'manual' renumbering 
    # is required as with the other schemes.
    
    # CDR1
    # CDR1 has a range from 27 (inc.) to 39 (exc.) and has a theoretical maximum length of 12.
    cdr1seq    = "".join([ x[1] for x in _regions[1] if x[1] != "-" ])
    cdr1length = len(cdr1seq) 
    si = 0
    prev_state = 26
    for ann in get_imgt_cdr(cdr1length, 12, 27, 39):
        if not ann:
            _numbering[1].append( ((prev_state+1, ' '), '-') )
            prev_state += 1
        else:
            _numbering[1].append( (ann, cdr1seq[si]) )
            prev_state = ann[0]
            si += 1

    # CDR2 
    # CDR2 has a range from 56 (inc.) to 66 (exc.) and has a theoretical length of 10.
    cdr2seq    = "".join([ x[1] for x in _regions[3] if x[1] != "-" ])
    cdr2length = len(cdr2seq)
    si = 0
    prev_state = 55
    for ann in get_imgt_cdr(cdr2length, 10, 56, 66):
        if not ann:
            _numbering[3].append( ((prev_state+1, ' '), '-') )
            prev_state += 1
        else:
            _numbering[3].append( (ann, cdr2seq[si]) )
            prev_state = ann[0]
            si += 1

    # FW3. We allow the HMM to place insertions. Technically all insertion points are taken care of but in reality insertions can
    # and do occur. No specification of where the insertions should be placed.


    # CDR3 
    # CDR3 has a range from 105 (inc.) to 118 (exc.). Insertions are placed on 112 and 111 symetrically. IMGT has a technical 
    # maximum length of 65 (13 positions, 26*2 insertions) . In practice ANARCI will not recognise CDR3s of this length.
    cdr3seq    = "".join([ x[1] for x in _regions[5] if x[1] != "-" ])
    cdr3length = len(cdr3seq)
    if cdr3length > 117: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    si = 0
    previous_state_id = 104
    for ann in get_imgt_cdr(cdr3length, 13, 105, 118):
        if ann is None:
            _numbering[5].append( ((previous_state_id+1, " "), "-"   ) )
            previous_state_id+=1
        else:
            _numbering[5].append( (ann, cdr3seq[si] ) )
            previous_state_id = ann[0]
            si+=1
  
    # Return the full vector and the start and end indices of the numbered region of the sequence
    return gap_missing( _numbering ), startindex, endindex
  
def get_imgt_cdr(length, maxlength, start, end):
    """
    Symmetrically number a CDR loop (e.g. CDRL1/CDRH2 for IMGT)
    @param length:      Define the length of target CDR
    @param maxlength:   Define the theoretical limit (e.g. L1 = 12 for the IMGT scheme)
    @param start, end:  Start and end position numbers
    """
    annotations = [ None for _ in range(max(length, maxlength)) ]
    if length == 0:
        return annotations
    elif length == 1:
        annotations[0] = (start, ' ')
        return annotations

    front, back = 0, -1
    #az = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" 
    #za = "ZYXWVUTSRQPONMLKJIHGFEDCBA"

    az = alphabet[:-1]
    za = az[::-1]

    for i in range(min(length, maxlength)):
        if i % 2:
            annotations[back] = (end + back, " ")
            back -= 1
        else:
            annotations[front] = (start + front, " ")
            front += 1

    # Add insertions around the centre point
    centrepoint = [ i for i,v in enumerate(annotations) if v == None ]
    if not centrepoint:
        return annotations

    centre_left  = annotations[min(centrepoint)-1][0] # Get the index right before the first None
    centre_right = annotations[max(centrepoint)+1][0] # Get the index right after  the first None

    # For cases with an even max length
    if not maxlength % 2:
        frontfactor, backfactor = maxlength//2, maxlength//2
    # For cases with an odd max length
    else:
        frontfactor, backfactor = (maxlength//2)+1, maxlength//2

    for i in range(max(0, length-maxlength)):
        if not i % 2:
            annotations[back] = (centre_right, za[back + backfactor])
            back -= 1
        else:
            annotations[front] = (centre_left, az[front - frontfactor])
            front += 1
    
    return annotations


#######
# Aho # 
#######
# Heuristic regapping based on the AHo specification as detailed on AAAAA website. Gap order depends on the chain type
def number_aho(state_vector, sequence, chain_type):
    """    
    Apply the Aho numbering scheme
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in IMGT scheme, I is an insertion. (All X's for IMGT)

    XXXXXXX XXX XXXXXXXXXXXXXX XXXXXXXXXXXXXXXX XXXXXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXX
    AAAAAAA BBB CCCCCCCCCCCCCC DDDDDDDDDDDDDDDD EEEEEEEEEEEEEEE FFFFFFFFFFFFFFFFFFFF HHHHHHHHHHHHHHHH IIIIIIIIIIIII JJJJJJJJJJJJJ KKKKKKKKKKK


    Regions - (N.B These do not match up with any particular definition of CDR)
    A. EMPTY (now included in B)
    B. 1-10 inclusive. Indel occurs at 8
    C. 11-24 inclusive.
    D. 25-42 inclusive (deletion surround 28) 32-42 inclusive (deletions surround 36)
    E. 43-57 inclusive 
    F. 58-77 inclusive (deletions surround 63). Alpha chains have deletions at 74,75
    G. EMPTY (now included in H)
    H. 78-93 inclusive  gaps on 86 then 85, insertions on 85 linearly
    I. 94-106 inclusive
    J. 107-138 inclusive gaps on 123 symetrically.
    K. 139-149 inclusive.
    
    """
    
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string =  'BBBBBBBBBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFHHHHHHHHHHHHHHHHIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK'
#                     1         2             3               4              5                   7               8            9            10


    region_index_dict = dict( list(zip( "ABCDEFGHIJK", list(range(11)) )) )
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:0,
                         2:0,
                         3:0,
                         4:2,
                         5:2,
                         6:2,
                         7:2,
                         8:2,
                         9:2,
                         10:21}

    n_regions = 11
    
    exclude_deletions = [1,3,4,5,7,9]    
    
    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############

    _numbering = [ _regions[0], _regions[1], _regions[2],[], _regions[4], [], _regions[6], [], _regions[8],_regions[9],_regions[10] ]

    ##################################
    # Move the indel in fw 1 onto 8  #
    ##################################

    # Place indels on 8 
    # Find the first recognised residue and change the expected length of the stretch given the starting point. 
    # This prevents n terminal deletions being placed at 8 incorrectly.
    length = len( _regions[1] )
    if length > 0:
        start = _regions[1][0][0][0] 
        stretch_len = 10 - (start -1)
        if length > stretch_len: # Insertions are present. Place on 8
            annotations = [ (_," ") for _ in range(start,9) ] + [ (8,alphabet[_]) for _ in range( length - stretch_len ) ] + [(9," "),(10," ")]
        else:
            ordered_deletions = [(8," ")] + [(_," ") for _ in range(start, 11) if _ != 8]
            annotations = sorted( ordered_deletions[max(stretch_len-length, 0):] )
        _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in range(length) ]

    #########
    # CDR 1 # - divided in two parts in the Aho scheme.
    ######### - gaps at 28 depending on the chain type.

    # "VH domains, as well as the majority of the VA domains, have a one-residue gap in position 28, VK and VB domains a two-residue
    # gap in position 27 and 28."

    # We use the link below as the reference for the scheme.
    # https://www.bioc.uzh.ch/plueckthun/antibody/Numbering/Alignment.html

    # Some of the header lines in these images are offset by one (VH)! The gaps really are centered at 28 and 36
    # https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VK.html
    # https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VL.html
    # https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VH.html
    # https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VA.html
    # https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VB.html
    # https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VG.html
    # https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VD.html

    # We gap the CDR1 in a heuristic way using the gaps. 
    # This means that CDR1 gapping will not always be correct. For example if one grafts a Kappa CDR1 loop onto a Lambda framework
    # the gapping patter might now be incorrect. 
    # Not a fan of being so prescriptive.

    # The CDR1 region included here ranges from AHo 25 to AHo 42 inclusive
    
    # The order in which the two loops are gapped is dependent on the chain type (see alignments in URLs above).
    # Not all lengths are defined as not all lengths were crystallised in 2001 (or today). Where no example of the length was 
    # available the rule followed is to continue gapping the C terminal 'loop', then the N terminal 'loop', then 31 then the fw. 
    # In all cases I have commented where the gapping is undefined. Note that for alpha chains the gapping rules are inconsistent.

    _L = 28,36,35,37,34,38,27,29,33,39,32,40,26,30,25,31,41,42
    #                           |-> undefined by AHo. Gapping C terminal loop then N terminal then 31, then fw.                              
    _K = 28,27,36,35,37,34,38,33,39,32,40,29,26,30,25,31,41,42
    #                                 |-> undefined by AHo. Gapping C terminal loop then N terminal then fw.                              
    _H = 28,36,35,37,34,38,27,33,39,32,40,29,26,30,25,31,41,42 
    #                        |-> undefined by AHo. Gapping C terminal loop then N terminal then fw.                              
    #                            N.B. The header on the alignment image for PDB_VH is offset by 1!
    _A = 28,36,35,37,34,38,33,39,27,32,40,29,26,30,25,31,41,42    
    #                              |-> undefined by AHo. Gapping C terminal loop then N terminal then fw.     
    #                            N.B The gapping is inconsistent for alpha chains. I follow the paper's statement that most VA have
    #                                one gap at 28 and remove 28 and 27 before removing 40. 
    _B = 28,36,35,37,34,38,33,39,27,32,40,29,26,30,25,31,41,42
    #                              |-> undefined by AHo. Gapping C terminal loop then N terminal then 31, then fw.
    _D = 28,36,35,37,34,38,27,33,39,32,40,29,26,30,25,31,41,42
    #                         |-> undefined by AHo. Gapping C terminal loop then N terminal then 31, then fw.
    #                         N.B only two sequence patterns available.                                
    _G = 28,36,35,37,34,38,27,33,39,32,40,29,26,30,25,31,41,42
    #                         |-> undefined by AHo. Gapping C terminal loop then N terminal then 31, then fw.
    #                         N.B only one sequence patterns available. Delta copied.
    
    ordered_deletions = { 'L':_L,'K':_K, 'H':_H, 'A':_A, 'B':_B, 'D':_D, 'G':_G } 

    length = len( _regions[3] )

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[chain_type][ max(18-length, 0): ] ) ]

    # Insertions are not described in the AHo scheme but must be included as there is a significant number of CDRH1s that are 
    # longer than the number of positions.
    insertions = max( length-18 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    elif insertions > 0:
        # They are placed on residue 36 alphabetically.
        insertat = annotations.index( (36, ' ') )+1 # Always 12 
        assert insertat == 12, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (36, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in range(length) ]

    #########
    # CDR 2 #
    #########
    # Gaps are placed symetically at 63.
    # For VA a second gap is placed at 74 and 75 according to the text in the paper. However, all the reference sequences show a
    # gap at 73 and 74 see:
    #      https://www.bioc.uzh.ch/plueckthun/antibody/Sequences/Rearranged/PDB_VA.html 
    # and 
    #      https://www.bioc.uzh.ch/plueckthun/antibody/Numbering/Alignment.html
    # Either I am mis-interpreting the text in the paper or there is something a little inconsistent here...
    # Given that *all* the numbered examples show the VA gap at 73 and 74 on the AAAAA website I have decided to implement this. 
    #
    
    # This region describes 58 to 77 inclusive
    
    if chain_type == 'A':
        ordered_deletions = [74,73,63,62,64,61,65,60,66,59,67,58,68,69,70,71,72,75,76,77]
    else:
        ordered_deletions = [63,62,64,61,65,60,66,59,67,58,68,69,70,71,72,73,74,75,76,77]

    length = len(_regions[5])

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[ max(20-length, 0): ] ) ]

    # Insertions are not described in the AHo scheme but must be included.
    insertions = max( length-20 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    elif insertions > 0:
        # They are placed on residue 63 alphabetically.
        insertat = annotations.index( (63, ' ') )+1 # Always 6
        assert insertat == 6, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (63, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[5] = [ (annotations[i], _regions[5][i][1]) for i in range(length) ]

    #########
    # FW3   ############################################
    # Move deletions onto 86 then 85. Insertions on 85 #
    ####################################################
    ordered_deletions = [86,85,87,84,88,83,89,82,90,81,91,80,92,79,93,78]
    length=len( _regions[7] )

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[ max(16-length, 0): ] ) ]

    # Insertions are not described in the AHo scheme but must be included.
    insertions = max( length-16 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    elif insertions > 0:
        # They are placed on residue 85 alphabetically.
        insertat = annotations.index( (85, ' ') )+1 # Always 8
        assert insertat == 8, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (85, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[7] = [ (annotations[i], _regions[7][i][1]) for i in range(length) ]


    #########
    # CDR 3 #
    #########
    # Deletions on 123. 
    # Point of the Aho scheme is that they have accounted for all possible positions.
    # Assumption is that no more insertions will occur.... 
    # We'll put insertions on 123 linearly.(i.e.ABCDEF...) if they ever do.

    ordered_deletions = [123,124,122,125,121,126,120,127,119,128,118,129,117,130,116,131,115,132,114,133,113,134,112,135,111,
                         136,110,137,109,138,108,107]
    
    length=len( _regions[9] )

    annotations = [ (i, ' ') for i in sorted( ordered_deletions[ max(32-length, 0): ] ) ]

    # Insertions are not described in the AHo scheme but must be included.
    insertions = max( length-32 , 0 ) 
    if insertions > 26: 
        return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    elif insertions > 0:
        # They are placed on residue 123 alphabetically.
        insertat = annotations.index( (123, ' ') )+1 # Always 17
        assert insertat == 17, 'AHo numbering failed'  
        annotations = annotations[:insertat] + [ (123, alphabet[a]) for a in range( insertions ) ] + annotations[insertat:]

    _numbering[9] = [ (annotations[i], _regions[9][i][1]) for i in range(length) ]

    # AHo includes one extra position than IMGT in what it considers the variable domain for light chains. 
    #If the last state is 148 and there is at least one more residue left, then add the residue to the numbering. 
    numbering = gap_missing( _numbering )
    if len(numbering) > 0:
        if numbering[-1][0] == (148, ' ') and numbering[-1][1] != '-' and endindex+1 < len(sequence):
            numbering.append( ( (149, ' '), sequence[endindex+1]) )
            endindex +=1

    return numbering, startindex, endindex


###########
# Chothia #
###########

# Heavy chains
def number_chothia_heavy(state_vector, sequence):
    """
    Apply the Chothia numbering scheme for heavy chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Chothia scheme, I is an insertion.

    XXXXXXXXXI XXXXXXXXXXXXX XXXXXXXIIIIXX XXXXXXXXXXXXXXXXXX XXXIXIIXXXX XXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXX XXXXXXXXIIIXX XXXXXXXXXXX'
    1111111111 2222222222222 3333333333333 444444444444444444 55555555555 666666666666666666666666666666666666666 7777777777777 88888888888'     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Put the insertions at Chothia position 6
     2  -  Simple mapping (treat "I" states as inserts and not own match states)
     3  -  CDRH1 - 30 (inc) to 34 (exc) put insertions on 31
     4  -  Simple mapping (treat "I" states as inserts and not own match states)
     5  -  CDRH2 - 52 (inc) 58 (exc) put insertions on 52 
     6  -  Simple mapping (treat "I" states as inserts and not own match states)
     7  -  CDRH3 93 (inc) to 103 (exc) put insertion on 100
     8  -  Simple mapping (treat "I" states as inserts and not own match states)


    Regions 1,3,5 and 7 are renumbered
    
    """

    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [0,2,4,6] # Don't put deletions in these regions

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    
    ###############
    # Renumbering #
    ###############

    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]
    
    # Chothia H region 1 (index 0)
    # Insertions are placed at Chothia position 6.
    # Count how many we recognised as insertion by the hmm
    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 
    # We will place all insertion in this region at Chothia position 6.
    if insertions:
        start = _regions[0][0][0][0] # The starting Chothia number as found by the HMM (could easily start from 2 for example)
        # I have a feeling this may be a source of a bug in very unusual cases. Can't break for now. Will catch mistakes in a validate function. 
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in range(start, 7) ] + [ (6, alphabet[_]) for _ in range(insertions) ] + [(7," "),(8," "),(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in range(length) ]
    else:
        _numbering[0] = _regions[0]

    
    # CDR1 
    # Chothia H region 3 (index 2)
    # put insertions onto 31
    length = len( _regions[2] )
    insertions = max(length - 11, 0) # Pulled back to the cysteine as heavily engineered cdr1's are not playing nicely

    if insertions:
        annotations = [(_, " ") for _ in range(23,32)] + [(31, alphabet[i]) for i in range(insertions) ] + [(32," "),(33," ")]
    else:
        annotations = [(_, " ") for _ in range(23,32)][:length-2] + [(32," "),(33," ")][:length]

    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in range(length) ]
 
    # CDR2
    # Chothia H region 5 (index 4) 
    # put insertions onto 52
    length = len( _regions[4] )
    # 50 to 57 inclusive
    insertions = max(length - 8, 0) # Eight positions can be accounted for, the remainder are insertions
    # Delete in the order, 52, 51, 50,53, 54 ,55, 56, 57
    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in range(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]

    # FW3 - insertions are annotated on 82. The first three are normal positions and annotated automatically. 
    # Additional insertions do not occur with the kabat or the chothia numbering scheme.
    # It does not make sense to place more than A, B, C on 82 as Martin and AHo work show that this is not a place that accepts 
    # additional insertions. 
    # The decision here is to allow the alignment to place additional insertions. This is in contrast to Martin where the region
    # is renumbered to place insertions on 72. 
     
    # CDR3
    # Chothia H region 7 (index 6) 
    # put insertions onto 100
    length = len( _regions[6] )    
    if length > 36: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="heavy")
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in range(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return gap_missing( _numbering ), startindex, endindex                                     

# Light chains
def number_chothia_light(state_vector, sequence):
    """
    Apply the Chothia numbering scheme for light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Chothia scheme, I is an insertion.
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX IIIIIIX XXXXXXXXXXXXXXXXXXXX XIIIIIIIXXX XXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXX XXXXXIIIIXX XXXXXXXXXXXXX
    11111111111111111111111111111 2222222 33333333333333333333 44444444444 5555555555555555555555555555555555555 66666666666 7777777777777
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 24 (inc) to 35 (exc) put insertions on 30
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 51 (inc) 55 (exc) put insertions on 52
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRL3 89 (inc) to 98 (exc) put insertion on 95
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    Region 2, 3 and 5 are renumbered
    
    """
    
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1: 0,
                         2:-6,
                         3:-6,
                         4:-13,
                         5:-16,
                         6:-20,
                         }    

    
    n_regions = 7
    
    exclude_deletions = [1,3,4,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    _numbering = [ _regions[0], [], _regions[2], [], _regions[4], [], _regions[6] ]
    
    
    ###############
    # Renumbering #
    ###############

    # CDR1 
    # Chothia L region 2 (index 1)
    # put insertions onto 30
    length = len( _regions[1] )
    insertions = max(length - 11, 0) # Eleven positions can be accounted for, the remainder are insertions
    # Delete forward from 31 
    annotations  =  [(24, " "),(25, " "), (26, " "), (27, " "), (28, " "),(29, " "),(30, " ")][:max(0,length)] 
    annotations += [(30, alphabet[i]) for i in range(insertions) ]
    annotations += [(31, " "),(32, " "),(33, " "),(34, " ")][ abs( min(0,length-11) ):] 
    _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in range(length) ]


    # CDR2
    # Chothia L region 4 (index 3) 
    # put insertions onto 52. 
    length = len( _regions[3] )
    insertions = max( length - 4, 0 )
    if insertions > 0:
        annotations  = [(51, " "),(52, " ")] + [(52, alphabet[i]) for i in range(insertions) ] + [(53, " "),(54, " ")]
        _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in range(length) ]
    else: # How to gap L2 in Chothia/Kabat/Martin is unclear so we let the alignment do it.
        _numbering[3] = _regions[3]
    
    # FW3
    # Insertions on 68. First deletion 68. Otherwise default to alignment
    length = len( _regions[4] )
    insertions = max(length - 34, 0)
    if insertions > 0: # Insertions on 68
        annotations = [(i," ") for i in range(55,69)]+[(68, alphabet[i]) for i in range(insertions) ]+[(i," ") for i in range(69,89)]
        _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]
    elif length == 33: # First deletion on 68
        annotations = [(i," ") for i in range(55,68)]+[(i," ") for i in range(69,89)]            
        _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]
    else: # More deletions - allow alignment to place them 
        _numbering[4] = _regions[4]


    # CDR3
    # Chothia L region 6 (index 5) 
    # put insertions onto 95
    length = len( _regions[5] )    

    if length > 35: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="light")
    _numbering[5]  = [ (annotations[i], _regions[5][i][1]) for i in range(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence

    return gap_missing( _numbering ), startindex, endindex    


#########
# Kabat #
#########

# Heavy chains
def number_kabat_heavy(state_vector, sequence):
    """
    Apply the Kabat numbering scheme for heavy chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Kabat scheme, I is an insertion.
    XXXXXXXXXI XXXXXXXXXXXXXXXXXXXX IIIIXXXXXX XXXXXXXXXXXXXXXX XIXII XXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXX XXXXXXIII XXXXXXXXXXXXX
    1111111111 22222222222222222222 3333333333 4444444444444444 55555 666666666666666666666666666666666666666666666 777777777 8888888888888
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Put the insertions at Chothia position 6
     2  -  Simple mapping (treat "I" states as inserts and not own match states)
     3  -  CDRH1 - 30 (inc) to 36 (exc) put insertions on 35
     4  -  Simple mapping (treat "I" states as inserts and not own match states)
     5  -  CDRH2 - 52 (inc) 58 (exc) put insertions on 52 
     6  -  Simple mapping (treat "I" states as inserts and not own match states)
     7  -  CDRH3 93 (inc) to 103 (exc) put insertion on 100
     8  -  Simple mapping (treat "I" states as inserts and not own match states)

    """
 
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [2,4,6]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 0, 2, 4, 6 regions in Chothia heavy
    
    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]


    # Kabat H region 1 (index 0)
    # Insertions are placed at Kabat position 6.
    # Count how many we recognised as insertion by the hmm
    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 
    # We will place all insertion in this region at Kabat position 6.
    if insertions:
        start = _regions[0][0][0][0] # The starting Kabat number as found by the HMM (could easily start from 2 for example)
        # I have a feeling this may be a source of a bug in very unusual cases. Can't break for now. Will catch mistakes in a validate function. 
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in range(start, 7) ] + [ (6, alphabet[_]) for _ in range(insertions) ] + [(7," "),(8," "),(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in range(length) ]
    else:
        _numbering[0] = _regions[0]
    
    
    # CDR1 
    # Kabat H region 3 (index 2)
    # Put insertions onto 35. Delete from 35 backwards
    length = len( _regions[2] )
    insertions = max(0,length - 13)
    annotations = [(_,' ') for _ in range(23, 36)][:length] 
    annotations += [(35, alphabet[i]) for i in range(insertions) ]
    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in range(length) ]
 
    # CDR2
    # Chothia H region 5 (index 4) 
    # put insertions onto 52
    length = len( _regions[4] )
    # 50 to 57 inclusive
    insertions = max(length - 8, 0) # Eight positions can be accounted for, the remainder are insertions
    # Delete in the order, 52, 51, 50,53, 54 ,55, 56, 57
    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in range(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]

    # FW3 - insertions are annotated on 82. The first three are normal positions and annotated automatically. 
    # Additional insertions do not occur with the kabat or the chothia numbering scheme.
    # It does not make sense to place more than A, B, C on 82 as Martin and AHo work show that this is not a place that accepts 
    # additional insertions. 
    # The decision here is to allow the alignment to place additional insertions. This is in contrast to Martin where the region
    # is renumbered to place insertions on 72. 
     
    # CDR3
    # Chothia H region 7 (index 6) 
    # put insertions onto 100
    length = len( _regions[6] )    
    if length > 36: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="kabat", chain_type="heavy") #  Chothia and Kabat the same here
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in range(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return gap_missing( _numbering ), startindex, endindex            
           
# Light chains    
def number_kabat_light(state_vector, sequence):
    """
    Apply the Kabat numbering scheme for light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Kabat scheme, I is an insertion.
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX IIIIIIX XXXXXXXXXXXXXXXXXXXX XIIIIIIIXXX XXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXX XXXXXIIIIXX XXXXXXXXXXXXX
    11111111111111111111111111111 2222222 33333333333333333333 44444444444 5555555555555555555555555555555555555 66666666666 7777777777777
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 24 (inc) to 35 (exc) put insertions on 27
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 51 (inc) 55 (exc) put insertions on 52
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRL3 89 (inc) to 96 (exc) put insertion on 95
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    """
    
    # Set up the numbering 
 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777'
    
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1: 0,
                         2:-6,
                         3:-6,
                         4:-13,
                         5:-16,
                         6:-20,
                         }    
    
    n_regions = 7
    
    exclude_deletions = [1,3,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    _numbering = [ _regions[0], [], _regions[2], [], _regions[4], [], _regions[6] ]
    
    
    ###############
    # Renumbering #
    ###############
    
    # CDR1 
    # Kabat L region 2 (index 1)
    # put insertions onto 27
    length = len( _regions[1] )
    insertions = max(length - 11, 0) # Eleven positions can be accounted for, the remainder are insertions
    # Delete forward from 28 
    annotations  =  [(24, " "),(25, " "), (26, " "), (27, " ")][:max(0,length)] 
    annotations += [(27, alphabet[i]) for i in range(insertions) ]
    annotations += [(28, " "),(29, " "),(30, " "),(31, " "),(32, " "),(33, " "),(34, " ")][ abs( min(0,length-11) ):] 
    _numbering[1] = [ (annotations[i], _regions[1][i][1]) for i in range(length) ]
  
    # CDR2
    # Chothia L region 4 (index 3) 
    # put insertions onto 52. 
    length = len( _regions[3] )
    insertions = max( length - 4, 0 )
    if insertions > 0:
        annotations  = [(51, " "),(52, " ")] + [(52, alphabet[i]) for i in range(insertions) ] + [(53, " "),(54, " ")]
        _numbering[3] = [ (annotations[i], _regions[3][i][1]) for i in range(length) ]
    else: # How to gap L2 in Chothia/Kabat/Martin is unclear so we let the alignment do it.
        _numbering[3] = _regions[3]


    # FW3
    # All insertions are placed by alignment. This is in contrast to Martin (and Chothia) where they are placed on 68. 
    # The kabat scheme was defined using a sequence alignment alone. In keeping with this, insertions in FW3 are also only placed
    # with respect to the sequence alignment (the HMM).

    # CDR3
    # Chothia L region 6 (index 5) 
    # put insertions onto 95
    length = len( _regions[5] )    

    if length > 35: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="kabat", chain_type="light")
    _numbering[5]  = [ (annotations[i], _regions[5][i][1]) for i in range(length)  ]

    return gap_missing( _numbering ), startindex, endindex    




#############################
# Martin (extended Chothia) #
#############################

# Heavy chains
def number_martin_heavy(state_vector, sequence):
    """
    Apply the Martin (extended Chothia) numbering scheme for heavy chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Martin scheme, I is an insertion.
    XXXXXXXXXI XXXXXXXXXXXXXXXXXXXX IIIIXX XXXXXXXXXXXXXXXXXXXX XIXII XXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXX XXXXXXIII XXXXXXXXXXXXX
    1111111111 22222222222222222222 333333 44444444444444444444 55555 666666666666666666666666666666666666666666666 777777777 8888888888888
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Put the insertions at Chothia position 8
     2  -  Simple mapping (treat "I" states as inserts and not own match states)
     3  -  CDRH1 - 30 (inc) to 34 (exc) put insertions on 31
     4  -  Simple mapping (treat "I" states as inserts and not own match states)
     5  -  CDRH2 - 52 (inc) 58 (exc) put insertions on 52 
     6  -  Simple mapping (treat "I" states as inserts and not own match states)
     7  -  CDRH3 93 (inc) to 103 (exc) put insertion on 100
     8  -  Simple mapping (treat "I" states as inserts and not own match states)


    Regions 1,3,5 and 7 are renumbered
    
    """
 
    # Set up the numbering 
 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888'
    
    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6,"8":7}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:0, 
                         1:-1,
                         2:-1,
                         3:-5,
                         4:-5,
                         5:-8,
                         6:-12,
                         7:-15}    
    
    n_regions = 8
    
    exclude_deletions = [2,4,5,6]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    

    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 0, 2, 4, 6 regions in Chothia heavy
    
    _numbering = [ [], _regions[1] , [], _regions[3] , [], _regions[5], [], _regions[7] ]
    
    # Chothia H region 1 (index 0)
    # Insertions are placed at Chothia position 8.
    # Count how many we recognised as insertion by the hmm
    insertions = len( [ 1 for _ in _regions[0] if _[0][1] != " " ] ) 
    # We will place all insertion in this region at Chothia position 8.
    if insertions:
        start = _regions[0][0][0][0] # The starting Chothia number as found by the HMM (could easily start from 2 for example)
        # I have a feeling this may be a source of a bug in very unusual cases. Can't break for now. Will catch mistakes in a validate function. 
        length = len( _regions[0] )
        annotations = [ (_, " ") for _ in range(start, 9) ] + [ (8, alphabet[_]) for _ in range(insertions) ] + [(9," ")]
        _numbering[0] =  [ (annotations[i], _regions[0][i][1]) for i in range(length) ]
    else:
        _numbering[0] = _regions[0]

    
    # CDR1 
    # Chothia H region 3 (index 2)
    # put insertions onto 31
    length = len( _regions[2] )
    insertions = max(length - 11, 0) # Pulled back to the cysteine as heavily engineered cdr1's are not playing nicely
    if insertions:
        annotations = [(_, " ") for _ in range(23,32)] + [(31, alphabet[i]) for i in range(insertions) ] + [(32," "),(33," ")]
    else:
        annotations = [(_, " ") for _ in range(23,32)][:length-2] + [(32," "),(33," ")][:length]
    _numbering[2] = [ (annotations[i], _regions[2][i][1]) for i in range(length) ]
 
    # CDR2
    # Chothia H region 5 (index 4) 
    # put insertions onto 52
    length = len( _regions[4] )
    # 50 to 57 inclusive
    insertions = max(length - 8, 0) # Eight positions can be accounted for, the remainder are insertions
    # Delete in the order, 52, 51, 50,53, 54 ,55, 56, 57
    annotations  =  [(50, " "),(51, " "), (52, " ")][:max(0,length-5)]
    annotations += [(52, alphabet[i]) for i in range(insertions) ]
    annotations += [(53, " "),(54, " "),(55, " "),(56, " "),(57, " ")][ abs( min(0,length-5) ):]
    _numbering[4] = [ (annotations[i], _regions[4][i][1]) for i in range(length) ]

    # FW3 
    # Place all insertions on 72 explicitly.
    # This is in contrast to Chothia implementation where 3 insertions are on 82 and then further insertions are placed by the 
    # alignment
    # Gaps are placed according to the alignment. 
    length = len( _regions[5] )
    insertions = max(length - 35, 0)
    if insertions > 0: # Insertions on 72
        annotations = [(i,' ') for i in range(58,73)]+[(72, alphabet[i]) for i in range(insertions) ]+[(i,' ') for i in range(73,93)]
        _numbering[5] = [ (annotations[i], _regions[5][i][1]) for i in range(length) ]
    else: # Deletions - all alignment to place them. 
        _numbering[4] = _regions[4]

     
    # CDR3
    # Chothia H region 7 (index 6) 
    # put insertions onto 100
    length = len( _regions[6] )    
    if length > 36: return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = get_cdr3_annotations(length, scheme="chothia", chain_type="heavy")
    _numbering[6]  = [ (annotations[i], _regions[6][i][1]) for i in range(length)  ]

    # Return the full vector and the start and end indices of the numbered region of the sequence
    return gap_missing( _numbering ), startindex, endindex                                    

# Light chains
def number_martin_light(state_vector, sequence):
    """
    Apply the Martin numbering scheme for light chains
    
    Rules should be implemented using two strings - the state string and the region string. 

    There are 128 states in the HMMs. Treat X as a direct match in Martin scheme, I is an insertion.
    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX IIIIIIX XXXXXXXXXXXXXXXXXXXX XIIIIIIIXXX XXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXX XXXXXIIIIXX XXXXXXXXXXXXX
    11111111111111111111111111111 2222222 33333333333333333333 44444444444 5555555555555555555555555555555555555 66666666666 7777777777777
     
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 30 (inc) to 31 (exc) put insertions on 30
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 51 (inc) 55 (exc) put insertions on 52 
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRL3 89 (inc) to 96 (exc) put insertion on 95
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    Region 2, 3 and 5 are renumbered
    
    """
    
    # The Martin and Chothia specification for light chains are very similar. Martin is more explicit in the location of indels 
    # but unlike the heavy chain these are additional instead of changes to the Chothia scheme. Thus, Chothia light is implemented
    # as martin light.
    return number_chothia_light(state_vector,sequence)


###########
# Wolfguy #
###########
# The Wolfguy numbering scheme is an in-house scheme used at Roche. It has been described publicly in the paper:
# Prediction of VH-VL domain orientation for antibody variable domain modeling. Bujotzek A. et al. Protein 2015 83(4) 681-95
#
# It is similar in gapping as IMGT and is defined only for heavy and light antibody chains. 
# Unlike other schemes the numbering denotes both the chain (heavy 101-499, light 501-799) and the region (less than -50 framework
# greater than -50 CDR). All CDRs of length less than 50 can be handled without the need for insertion codes. Numbering of the 
# framework behaves similarly to IMGT in that all positions are assumed to be accounted for. Framework insertions are placed by 
# the alignment. 
#
# Numbering of all CDRs is performed symmetrically with the exception of CDRL1. In this case the CDR is numbered according to a 
# pattern specific to the canonical class. This is recognised by length and by sequence similarity to a consensus sequence. If a 
# length has not been observed it is numbered symmetrically. 


def number_wolfguy_heavy(state_vector, sequence):
    """
    Apply the wolfguy numbering scheme for heavy chains 

    The scheme numbers the sequence using different segments so that the numbering tells you
    where in the antibody the sequence is describing. 

    XXXXXXXXXIXXXXXXXXXXXXXXXX XXXXXXXXXXXXXX XXXXXXXXXXXXXX XXXXXXXXXXXXXXXXXXIX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXX
    11111111111111111111111111 22222222222222 33333333333333 44444444444444444444 555555555555555555555555555555 6666666666666 77777777777'
    
    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRH1 - 155-199 (inc). Gap symmetrically about 175-176.
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRH2 - 251-299 (inc). Gap symmetrically about 271-272, then gap back from 294.
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  CDRH3 331,332 and 351-399 (inc). Gap according to the  
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
    
     Start gaps on rhs each time.
    """
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '11111111111111111111111111222222222222223333333333333344444444444444444444555555555555555555555555555555666666666666677777777777'

    region_index_dict = {"1":0,"2":1,"3":2,"4":3,"5":4,"6":5,"7":6}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:100, 
                         1:124,
                         2:160,
                         3:196,
                         4:226,
                         5:244,
                         6:283}    
    
    n_regions = 7
    
    exclude_deletions = [1,3,5]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 1, 3, 5 regions in wolfguy heavy    
    _numbering = [ _regions[0], [] , _regions[2], [], _regions[4] , [], _regions[6] ]

    # CDRH1
    # Delete symmetrically about 177. Delete right first.
    # May have to change this to reflect where the point of symmetry is
    ordered_deletions = [151]
    for p1,p2 in zip( list(range(152,176)), list(range(199, 175,-1))): ordered_deletions += [ p1,p2 ]
    length = len( _regions[1] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[1]  = [ ((annotations[i]," "), _regions[1][i][1]) for i in range(length)  ]
    
    # CDRH2
    # Delete symmetrically about 271. Delete right first.
    # Then delete right from 288
    ordered_deletions = [251]
    for p1,p2 in zip( list(range(252,271)), list(range(290, 271,-1))): ordered_deletions += [ p1,p2 ]
    ordered_deletions.append( 271 )
    ordered_deletions = list(range( 299, 290, -1)) + ordered_deletions
    length = len( _regions[3] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[3]  = [ ((annotations[i]," "), _regions[3][i][1]) for i in range(length)  ]

    # CDRH3	
    # Delete symmetrically about 374. Delete right first. 
    # Scheme changes at length 8
    # Scheme changes at length 12
    ordered_deletions = []
    for p1,p2 in zip( list(range(356,374)), list(range(391, 373,-1))): ordered_deletions += [ p1,p2 ]
    ordered_deletions = [ 354, 394, 355, 393, 392 ] + ordered_deletions 
    ordered_deletions = [331,332] + [ 399, 398, 351, 352, 397, 353, 396, 395 ] + ordered_deletions
    length = len( _regions[5] )

    if length > len(ordered_deletions): return [], startindex, endindex # Too many insertions. Do not apply numbering.     
    annotations = sorted(ordered_deletions[:length])
    _numbering[5]  = [ ((annotations[i]," "), _regions[5][i][1]) for i in range(length)  ]    
  
    # Return the full vector and the start and end indices of the numbered region of the sequence
    return sum( _numbering, [] ), startindex, endindex   
            

def number_wolfguy_light(state_vector, sequence):
    """
    Apply the wolfguy numbering scheme for light chains 

    The scheme numbers the sequence using different segments so that the numbering tells you
    where in the antibody the sequence is describing. 

    XXXXXXX XXX XXXXXXXXXXXXX XXXXXXXXXXXXXXXXX XXXXXXXXXXXXXXX XXXXXXXXXXXXXX XXXIXXXXXXX XXXX XXXXXXXXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXX
    1111111 AAA BBBBBBBBBBBBB 22222222222222222 333333333333333 44444444444444 55555555555 6666 77777777777777777777 8888888888888 99999999999

    Regions - (N.B These do not match up with any particular definition of CDR)
     1  -  Simple mapping (treat "I" states as inserts and not own match states)
     A  -  Move indels onto 508
     B  -  Simple mapping (treat "I" states as inserts and not own match states)
     2  -  CDRL1 - 551-599 (inc). Assign via the matching consensus sequence and length.
     3  -  Simple mapping (treat "I" states as inserts and not own match states)
     4  -  CDRL2 - 651-699 (inc). Gap about 673 then right from 694
     5  -  Simple mapping (treat "I" states as inserts and not own match states)
     6  -  Move indels onto 713 and 714
     7  -  Simple mapping (treat "I" states as inserts and not own match states)
     8  -  CDRL3 751-799 (inc). Gap symmetrically about 374-375 
     9  -  Simple mapping (treat "I" states as inserts and not own match states)
    
    """
    # Set up the numbering 
 
    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
                    
    # Region string - regions that should be treated separately in putting the numbering together 
    region_string = '1111111AAABBBBBBBBBBBBB222222222222222223333333333333334444444444444455555555555666677777777777777777777888888888888899999999999'

    region_index_dict = {"1":0,"A":1,"B":2,"2":3,"3":4,"4":5,"5":6,"6":7,"7":8,"8":9,"9":10}
    
    # Define how the scheme's numbering differs from IMGT at the start of each region. 
    # This is updated in the loop below
    rels              =  {0:500,
                         1:500,
                         2:500,    
                         3:527,
                         4:560,
                         5:595,
                         6:631,
                         7:630,
                         8:630,                                                  
                         9:646,
                         10:683}    
    
    n_regions = 11
    
    exclude_deletions = [1,3,5,7,9]

    _regions, startindex, endindex = _number_regions(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)
    
    ###############
    # Renumbering #
    ###############
        
    # Renumbering required for 1, 3, 5 regions in wolfguy heavy    
    _numbering = [ _regions[0], [], _regions[2], [] , _regions[4], [], _regions[6], [], _regions[8], [], _regions[10] ]


    # Gaps in the first section go 508 instead of the imgt 510 equivalent
    length = len(_regions[1] )
    annotations = sorted([ (510,' '), (509, ' '), (508, ' ')][ :length ] + [(508,a) for a in alphabet[:max(0, length-3)]])  
    _numbering[1]  = [ (annotations[i], _regions[1][i][1]) for i in range(length)  ]
    
    # CDRL1
    # Number by predicting the canonical 
    length = len(_regions[3] )
    annotations = _get_wolfguy_L1( _regions[3], length)
    _numbering[3]  = [ ((annotations[i]," "), _regions[3][i][1]) for i in range(length)  ]
    
    # CDRL2
    # Delete about 673. Finally delete right from 694. Maintain 651 as the last deletion
    ordered_deletions = []
    for p1,p2 in zip( list(range(652,673)), list(range(694, 672,-1))): ordered_deletions += [ p2,p1 ]
    ordered_deletions = [651] + list(range( 699, 694, -1)) + ordered_deletions + [673]

    length = len( _regions[5] )
    annotations = sorted(ordered_deletions[:length])
    _numbering[5]  = [ ((annotations[i]," "), _regions[5][i][1]) for i in range(length)  ]


    # The placement of the indel in wolfguy is different to that in imgt
    length = len( _regions[7] )
    insertions = max( 0, length - 4 )
    annotations = [(711, ' '), (712, ' '), (713, ' '), (714, ' ')][:length] + [ (714, a) for a in alphabet[:insertions] ]    
    _numbering[7]  = [ (annotations[i], _regions[7][i][1]) for i in range(length)  ]  
    
    # CDRL3
    # Delete symmetrically about 775. Delete right first. Finally delete 798 and 799
    ordered_deletions = []
    for p1,p2 in zip( list(range(751,775)), list(range(799, 775,-1))): ordered_deletions += [ p1,p2 ]
    ordered_deletions.append( 775 )
  
    length = len( _regions[9] )
    if length > len(ordered_deletions): return [], startindex, endindex # Too many insertions. Do not apply numbering. 
    annotations = sorted(ordered_deletions[:length])
    _numbering[9]  = [ ((annotations[i]," "), _regions[9][i][1]) for i in range(length)  ]  
  
    # Return the full vector and the start and end indices of the numbered region of the sequence
    return sum( _numbering, [] ), startindex, endindex  


def _get_wolfguy_L1(seq, length):
    """
    Wolfguy's L1 annotation is based on recognising the length and the sequence pattern defined
    by a set of rules. If the length has not been characterised, we number symmetrically about the
    middle of the loop.
    """
    
    # These are the annotations for different lengths of L1 according to the wolfguy definitions.
    L1_sequences = {
    9: [['9',     'XXXXXXXXX', [551, 552, 554, 556, 563, 572, 597, 598, 599]]], 
    10: [['10',   'XXXXXXXXXX', [551, 552, 553, 556, 561, 562, 571, 597, 598, 599]]], 
    11: [['11a',  'RASQDISSYLA', [551, 552, 553, 556, 561, 562, 571, 596, 597, 598, 599]], 
         ['11b',  'GGNNIGSKSVH', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]], 
         ['11b.2','SGDQLPKKYAY', [551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599]]], 
    12: [['12a',  'TLSSQHSTYTIE', [551, 552, 553, 554, 555, 556, 561, 563, 572, 597, 598, 599]], 
         ['12b',  'TASSSVSSSYLH', [551, 552, 553, 556, 561, 562, 571, 595, 596, 597, 598, 599]], 
         ['12c',  'RASQSVxNNYLA', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]], 
         ['12d',  'rSShSIrSrrVh', [551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599]]], 
    13: [['13a',  'SGSSSNIGNNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 597, 598, 599]], 
         ['13b',  'TRSSGSLANYYVQ', [551, 552, 553, 554, 556, 561, 562, 563, 571, 572, 597, 598, 599]]], 
    14: [['14a',  'RSSTGAVTTSNYAN', [551, 552, 553, 554, 555, 561, 562, 563, 564, 571, 572, 597, 598, 599]], 
         ['14b',  'TGTSSDVGGYNYVS', [551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 596, 597, 598, 599]]], 
    15: [['15',   'XXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 594, 595, 596, 597, 598, 599]]], 
    16: [['16',   'XXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 594, 595, 596, 597, 598, 599]]], 
    17: [['17',   'XXXXXXXXXXXXXXXXX', [551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 584, 594, 595, 596, 597, 598, 599]]]
    }    

    if length in L1_sequences: # Use the pre-defined motif 
        # Find the maximum scoring canonical form for this length. 
        curr_max = None, -10000
        for canonical in L1_sequences[length]:
            sub_score = 0
            for i in range( length ):
                try:
                    sub_score += blosum62[ (seq[i][1].upper(), canonical[1][i].upper() ) ]
                except KeyError:
                    sub_score += blosum62[ (canonical[1][i].upper(), seq[i][1].upper() ) ]
            if sub_score > curr_max[1]:
                curr_max = canonical, sub_score

        # return the annotations
        return curr_max[0][2]
    else: # Use a symmetric numbering about the anchors.
        ordered_deletions = []
        for p1,p2 in zip( list(range(551,575)), list(range(599, 575,-1))): ordered_deletions += [ p2,p1 ]
        ordered_deletions.append(575)
        return sorted( ordered_deletions[:length] )

def gap_missing( numbering ):
    '''
    Place gaps when a number is missing. All except wolfguy are continuously numbered
    '''
    # Gaps placed where a number is not present
    num = [ ((0,' '),'-') ]
    for p, a in sum( numbering, [] ):  
        if p[0] > num[-1][0][0]+1:
            for _i in range( num[-1][0][0]+1, p[0] ):
                num.append( ((_i, ' '), '-' ) )
        num.append( (p,a) )
    return num[1:]


######################
# Annotation of CDR3 #
######################
    
def get_cdr3_annotations(length, scheme="imgt", chain_type=""):
    """
    Given a length of a cdr3 give back a list of the annotations that should be applied to the sequence.
    
    This function should be depreciated
    """
    az = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" 
    za = "ZYXWVUTSRQPONMLKJIHGFEDCBA"
    
    if scheme=="imgt":
        start, end = 105, 118 # start (inclusive) end (exclusive)
        annotations = [None for _ in range(max(length,13))]
        front = 0
        back  = -1
        assert (length-13) < 50, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        for i in range(min(length,13)):
            if i%2:
                annotations[back] = (end+back, " ")
                back -= 1
            else:
                annotations[front] = (start+front, " ")
                front += 1
        for i in range(max(0,length-13)): # add insertions onto 111 and 112 in turn
            if i%2:
                annotations[back] = (112, za[back+6])
                back-=1
            else:
                annotations[front] = (111, az[front-7])
                front +=1        
        return annotations

    elif scheme in [ "chothia", "kabat"] and chain_type=="heavy": # For chothia and kabat
        # Number forwards from 93
        insertions = max(length - 10, 0)
        assert insertions < 27, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        ordered_deletions = [ (100, ' '), (99,' '), (98,' '), (97,' '), (96,' '), (95,' '), (101,' '),(102,' '),(94,' '), (93,' ') ]
        annotations = sorted( ordered_deletions[ max(0, 10-length): ] + [ (100,a) for a in az[:insertions ] ] )
        return annotations

    elif scheme in [ "chothia", "kabat"] and chain_type=="light":
        # Number forwards from 89
        insertions = max(length - 9, 0)
        assert insertions < 27, "Too many insertions for numbering scheme to handle" # We ran out of letters.
        ordered_deletions = [ (95,' '),(94,' '),(93,' '),( 92,' '),(91,' '),(96,' '),(97,' '),(90,' '),(89,' ') ]
        annotations = sorted( ordered_deletions[ max(0, 9-length): ] + [ (95,a) for a in az[:insertions ] ] )
        return annotations

    else:
        raise AssertionError("Unimplemented scheme")

