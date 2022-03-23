"""
Program to reformat the v and j alignments for each species into a single alignment
with combinations of the two genes.

v genes will be trimmed to imgt position 107 ( three after the cysteine at 104 )

j genes are a little bit more complicated as they are not aligned.

We therefore align them with muscle and *do not allow* gaps to open. Use a high gapopen
penalty.

-8 seems to be the first integer to work. Set it to -10.
muscle -in js.fasta -gapopen -10

We have a reference sequence in the alignment for which the numbering is known. Therefore
make the right hand most column imgt 128. Any further are trimmed.

Once aligned, all residue before imgt 116 are ignored (101 in Chothia heavy).
All before are changed to gaps. The rationale is the alignment in this region is
not structurally valid anyway. Want the columns because of the potential imgt states. 

HMMER needs a single alignment to build a HMM. We therefore make putative germline sequences
by combining all the v sequences with all the j sequences. (Each v sequence appears n(j) times).

These are chucked into curated alignments and hmmbuild used to create the HMMs. 
   
"""

import os, sys
from subprocess import Popen, PIPE
from FastaIO import chunkify

amino_acids = sorted(list("QWERTYIPASDFGHKLCVNM"))
acid_set = set( amino_acids+["."])

file_path  = os.path.split(__file__)[0]
fasta_path = os.path.join( file_path, "IMGT_sequence_files", "fastafiles" )
curated_path = os.path.join( file_path, "curated_alignments" )

all_species = ["Homo_sapiens",
           "Mus",
           "Rattus_norvegicus",
           "Oryctolagus_cuniculus",
           "Macaca_mulatta",
           "Sus_scrofa",
           "Vicugna_pacos",
           "Bos_taurus"]

all_tr_species = ["Homo_sapiens",
           "Mus",
]
translations = {"Homo_sapiens":"human",
           "Mus":"mouse",
           "Rattus_norvegicus":"rat",
           "Oryctolagus_cuniculus":"rabbit",
           "Macaca_mulatta":"rhesus",
           "Sus_scrofa":"pig",
           "Vicugna_pacos":"alpaca",
           "Bos_taurus":"cow"}


def read_alignment(input_file, read_all=False, region_name=""):
    """
    """
    imgt_fields =  ["accession_number",
    "allele",  
    "species",  
    "functionality",  
    "region",  
    "start_and_end_positions_IMGT/LIGM-DB_accession_number", 
    "number_of_nucleotides", 
    "codon_start", 
    "number_nucleotides_added_in_5'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_added_in_3'_compared_IMGT/LIGM-DB", 
    "number_nucleotides_to_correct_sequencing_errors", 
    "number_of_amino_acids", 
    "number_of_characters", 
    "partial",  
    "reverse"]

    records = {}
    try:
        handle = open(input_file, "r")
    except IOError:
        print('Warning file', input_file, 'could not be found')
        return records

    region=""
    for record in chunkify(handle):

        fields = dict(list(zip( imgt_fields, record.description.split("|"))) )
        sequence = record.seq 
        # These are the ones we care about and will be used
        try:
            if fields['accession_number'] == 'None':continue
            if fields["functionality"]=="F" and not fields["partial"].strip() and not fields["reverse"].strip():
                if read_all:
                    pass
                elif fields["allele"].split("*")[-1].strip()!="01": 
                    continue
                if set(list(sequence))- acid_set:
    #                print >> sys.stderr,"Unexpected character in sequence"
    #                print >> sys.stderr,sequence
                    continue

                if fields["region"] == region_name:
                    records[ (fields["species"], fields[ "allele" ] ) ] = sequence
                elif region_name.startswith("C"):
                    if len(sequence) < 100: 
                        continue # Filter out partial sequences that IMGT have not....
                elif region:
                    assert fields["region"]==region, "The region for some the entries is different"

                region=fields["region"]
                records[ (fields["species"], fields[ "allele" ] ) ] = sequence
        except KeyError:
            print("Something wrong with the file %s"%input_file)
            continue
            
    handle.close()
    return records

def read_fasta(filename):
    """
    Read a sequence file and parse as description, string 
    """
    handle = open(filename, "r")
    records = [ [s.description, s.seq.replace(" ","")] for s in chunkify(handle) ]
    handle.close()
    return records

def write_fasta( sequences ):
    """
    Write a fasta file containing all sequences
    """   
    filename = os.path.join( file_path, "muscle_alignments", "all_js.fasta" )
    with open( filename, "w" ) as outfile:
        for al in sequences:
            for s in sequences[al]:
                print(">%s|%s|%s|%s"%tuple( list(al)+list(s) ), file=outfile)
                print(sequences[al][s], file=outfile)
    return filename


def format_c_genes(calignments, gene_name=""):    
    new_calignments = {} 
    for entry in calignments:
        if len(calignments[entry]) == 0:
            continue
        new_calignments[entry] = {}
        for seq in calignments[entry]:

            cspecies, callele = seq
            sequence = calignments[entry][seq]
                        
            # IMGT has a different alignment for C-genes than V and J genes.
            # This condition filters out the two sequeces that are not in the consistent format for C gene.
            #tttt = sequence[:132].ljust( 132 ).replace(" ",".")
            tttt = sequence[:149].ljust( 149 ).replace(" ",".")
            if tttt[104] != "C" or tttt[33] != "C": 
                print("Something wrong with ", entry, gene_name, sequence)
                continue

            max_length = 149
            new_name = "%s_%s_%s" % (gene_name, cspecies, callele)
            new_calignments[entry][ new_name ] = sequence[:max_length].ljust( max_length ).replace(" ",".")

    return new_calignments

def format_j_genes(jalignments):

    reference = ("WFAYWGQGTLVTVSA", 4  , 19 )
    #                 seq           start  end

    ffile = write_fasta(jalignments)
    al_filename = os.path.join( file_path, "muscle_alignments", "all_js_aligned.fasta" )
    
    if sys.platform == "darwin":
        pr = Popen( [ "muscle_macOS", "-in", ffile, "-gapopen", "-10", "-out", al_filename, ], stdout=PIPE, stderr=PIPE )
    else:
        pr = Popen( [ "muscle", "-in", ffile, "-gapopen", "-10", "-out", al_filename, ], stdout=PIPE, stderr=PIPE )
    o, e = pr.communicate()
    aligned = read_fasta( al_filename )
    new_jalignments = {} 

    # Find the reference sequence and what we need to do to map
    for name, sequence in aligned:
        if name == "Mus|H|Mus musculus|IGHJ3*01":
            ref_aligned = sequence
            break
    start = ref_aligned.index( reference[0] )
    if start > reference[1]:
        START = start+1-reference[1]
    else:
        START = 0
    END = start + 15


    for name, sequence in aligned:
        species, chain_type, id1, id2 =  name.strip(">").split("|")

        if (species, chain_type) not in new_jalignments:   
            new_jalignments[(species, chain_type)] = {}
        # We take the last 13 of the new alignment and pad into 20 long string 
        new_jalignments[(species, chain_type)][ (id1, id2) ] = sequence[START: END][-14:].rjust(20).replace(" ", ".")
    return new_jalignments

def format_v_genes(valignments):
    """
    Take upto and including imgt position 108 in the alignment. Pad with gaps on the right side
    """
    
    new_valignments = {} 
    for entry in valignments:
        species, chain_type = entry 
        new_valignments[entry] = {}
        for seq in valignments[entry]:
            sequence = valignments[entry][seq]
            if chain_type == "L" and translations[species] == "rhesus":
                sequence = rhesus_lambda(sequence)
            elif chain_type == "A" and translations[species] == "mouse":
                sequence = mouse_alpha(sequence)
            elif chain_type == "D" and translations[species] == "mouse":
                sequence = mouse_delta(sequence)
            new_valignments[entry][ seq ] = sequence[:108].ljust( 108 ).replace(" ",".")
            if new_valignments[entry][ seq ][103] != "C" or new_valignments[entry][ seq ][22] != "C": 
                sys.stderr.write("Warning - this alignment doesn't feature CYS at position 23 and/or position 104.\n")
                sys.stderr.write("%s,%s,%s\n" % (new_valignments[entry][ seq ], entry, seq))

    return new_valignments

def mouse_delta(sequence):
    """
    Mouse delta chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.

    This is particularly bad because alignment is not even consistent within the chain and species!!!

    Remove and return
    """
    # Check in here because not all are bad...recheck again in the format v genes just to make sure.
    if sequence[103] != "C" or sequence[22] != "C":
        return sequence[ : 8 ] + sequence[ 9:85 ] + sequence[86:]
    return sequence

def rhesus_lambda(sequence):
    """
    Rhesus lambda chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
    Remove and return
    """
    return sequence[:20]+sequence[21:51]+ sequence[53:] 

def mouse_alpha(sequence):
    """
    Mouse alpha chains have insertions. This screws up the alignment to everything else - not really IMGT gapped.
    Remove and return
    """
    return sequence[:8]+sequence[9:85]+sequence[86:]

def combine_sequences(vsequences, jsequences):
    """
    Do a pairwise combination of the v and j sequences to get putative germline sequences for the species.
    
    """
    combined_sequences = {}
    for v in vsequences:
        vspecies, vallele = v
        for j in jsequences:
            _, jallele= j
            combined_sequences[("%s_%s_%s"%(vspecies, vallele,jallele)).replace(" ", "_")] = vsequences[v] + jsequences[j]
    return combined_sequences       

def make_putative_alignments( vsequences, jsequences, calignments = None ):
    all_sequences = {}
    for species, chain_type in vsequences:
        if (species, chain_type) not in jsequences or (species, chain_type) not in vsequences: continue
        combined_sequences = combine_sequences( vsequences[ (species, chain_type) ], jsequences[ (species, chain_type) ] )
        all_sequences[ (species, chain_type) ] = combined_sequences
        output_stockholm( combined_sequences, "%s_%s"%(translations[species], chain_type) )

    # Write just the V and J combinations
    output_stockholm_all( all_sequences )
     
    # Write the V and J combinations and the c-domains
    if calignments:
        output_stockholm_all_and_C(all_sequences, calignments )


def write_germlines(vsequences, jsequences):
    """
    Compile a dictionary containing all the v and j germline sequences.
    """

    all_gene_alignments = {"J":{},"V":{}}

    for species, chain_type in vsequences:
        for ((_,gene), seq) in vsequences[ (species, chain_type) ].items():
            assert len(seq)==108, species+_+gene+chain_type+_+seq+str(len(seq))
            try:
                all_gene_alignments["V"][ chain_type ][ translations[species] ][ gene ] = seq.replace(".","-") + "-"*20
            except KeyError:        
                try:
                    all_gene_alignments["V"][ chain_type ][ translations[species] ]= { gene : seq.replace(".","-") + "-"*20 }
                except KeyError:
                    try:
                        all_gene_alignments["V"][ chain_type ] = { translations[species] : { gene : seq.replace(".","-") + "-"*20 } } 
                    except KeyError:
                        all_gene_alignments["V"] = { chain_type : { translations[species] : { gene : seq.replace(".","-") + "-"*20 } } }

        for ((_,gene), seq) in jsequences.get((species, chain_type),{} ).items():
            assert len(seq)==20
            try:
                all_gene_alignments["J"][ chain_type ][ translations[species] ][ gene ] = "-"*108 + seq.replace(".","-")
            except KeyError:        
                try:
                    all_gene_alignments["J"][ chain_type ][ translations[species] ]= { gene : "-"*108 + seq.replace(".","-") }
                except KeyError:
                    try:
                        all_gene_alignments["J"][ chain_type ] = { translations[species] : { gene : "-"*108 + seq.replace(".","-")} } 
                    except KeyError:
                        all_gene_alignments["J"] = { chain_type : { translations[species] : { gene : "-"*108 + seq.replace(".","-") } } }
        output_python_lookup(all_gene_alignments)


def output_python_lookup(all_gene_alignments, path=None):
    """
    Format a lookup table for the germline sequences. This can then be used by the final program.
    """

    if path is None:
        path = curated_path
    filename = os.path.join( path, "germlines.py")
    with open(filename,'w') as outfile:
        print("all_germlines = "+repr(all_gene_alignments), file=outfile)

def write_stockholm( sequences, ID, outfile):
        print("# STOCKHOLM 1.0", file=outfile)
        print("#=GF ID %s"%ID, file=outfile)
        
        pad_length = max(list(map(len, list(sequences.keys()))))+1
        for s in sequences:
            print(s.replace(" ", "_").ljust(pad_length), sequences[s].replace(".","-"), file=outfile)
        print("#=GC RF".ljust(pad_length), "x"*len(sequences[s]), file=outfile)
        print("//", file=outfile)


def output_C_alignments(alignments, c_name):
    """
    Write a stockholm for all C domains. 
    """
    for species, chain_type in alignments:
        output_stockholm( alignments[(species, chain_type) ], "%s_%s_%s"%(c_name, translations[species], chain_type) )

def output_stockholm_all_and_C(all_sequences, all_C_alignments, path=None):
    """
    Output a minimal stockholm alignment file for all sequences. 
    """
    if path is None:
        path = curated_path

    filename = os.path.join( path, "ALL_AND_C.stockholm")
    with open( filename, "w") as outfile:
        for species, chain_type in all_sequences:
            sequences = all_sequences[(species, chain_type)]
            l = len(list(sequences.values())[0])
            assert all( [1 if l == len(sequences[s]) else 0 for s in sequences]), "Not all sequences in alignment are the same length"
            write_stockholm( sequences, "%s_%s"%(translations[species], chain_type), outfile)

        for c_name in all_C_alignments:
            for species, chain_type in all_C_alignments[c_name] :
                write_stockholm( all_C_alignments[c_name][(species, chain_type) ], "%s_%s_%s"%(translations[species], chain_type,c_name), outfile)
    return filename      


def output_stockholm_all(all_sequences, path=None):
    """
    Output a minimal stockholm alignment file for all sequences. 
    """
    if path is None:
        path = curated_path

    filename = os.path.join( path, "ALL.stockholm")
    with open( filename, "w") as outfile:
        for species, chain_type in all_sequences:
            sequences = all_sequences[(species, chain_type)]
            l = len(list(sequences.values())[0])
            assert all( [1 if l == len(sequences[s]) else 0 for s in sequences]), "Not all sequences in alignment are the same length"
            write_stockholm( sequences, "%s_%s"%(translations[species], chain_type), outfile)

    return filename      

def output_stockholm(sequences, name, path=None):
    """
    Output a minimal stockholm alignment file. 
    """
    if path is None:
        path = curated_path

    filename = os.path.join( path, "%s.stockholm"%name)
    l = len(list(sequences.values())[0])
    
    assert all( [1 if l == len(sequences[s]) else 0 for s in sequences]), "Not all sequences in alignment are the same length"
    
    with open( filename, "w") as outfile:
        write_stockholm( sequences, name, outfile)

    
    return filename      

def main():
    """
    Read in the raw v and j alignments
    Format them and combine the sequences
    """
    print("\nFormatting alignments\n")
    valignments, jalignments = {},{}
    all_valignments, all_jalignments = {},{}
    #ccalignments, c1alignments, c2alignments, c3alignments = {}, {}, {}, {}

    print("IGs")
    for species in all_species:
        for chain_type in "HKL":
            if not os.path.isfile( os.path.join( fasta_path, "%s_%sV.fasta" % (species,chain_type)) ):
                continue

            print(species, chain_type)
            valignments[ (species, chain_type) ]  = read_alignment( os.path.join( fasta_path , "%s_%sV.fasta"%(species, chain_type) ), region_name= "V-REGION" )
            jalignments[ (species, chain_type) ]  = read_alignment( os.path.join( fasta_path , "%s_%sJ.fasta"%(species, chain_type) ), region_name= "J-REGION") 

            ### Comment out if you want constant regions?
            #if chain_type == "H":
            #    c1alignments[ (species, chain_type) ]   = read_alignment( os.path.join( fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="CH1") 
            #    c2alignments[ (species, chain_type) ]   = read_alignment( os.path.join( fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="CH2") 
            #    c3alignments[ (species, chain_type) ]   = read_alignment( os.path.join( fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="CH3") 
            #else:
            #    ccalignments[ (species, chain_type) ]   = read_alignment( os.path.join( fasta_path , "%s_%sC.fasta"%(species, chain_type) ), region_name="C-REGION") 

            all_valignments[ (species, chain_type) ]  = read_alignment( os.path.join( fasta_path , "%s_%sV.fasta"%(species, chain_type) ), region_name= "V-REGION", read_all=True)
            all_jalignments[ (species, chain_type) ]  = read_alignment( os.path.join( fasta_path , "%s_%sJ.fasta"%(species, chain_type) ), region_name= "J-REGION", read_all=True) 

    print("\nTRs")
    for species in all_tr_species:
        for chain_type in "ABGD":
            if not os.path.isfile( os.path.join( fasta_path, "%s_%sV.fasta" % (species,chain_type)) ):
                continue

            print(species, chain_type)
            valignments[ (species, chain_type) ]   = read_alignment( os.path.join( fasta_path , "%s_%sV.fasta"%(species, chain_type) ))
            jalignments[ (species, chain_type) ]   = read_alignment( os.path.join( fasta_path , "%s_%sJ.fasta"%(species, chain_type) )) 
            all_valignments[ (species, chain_type) ]  = read_alignment( os.path.join( fasta_path , "%s_%sV.fasta"%(species, chain_type) ), read_all=True)
            all_jalignments[ (species, chain_type) ]  = read_alignment( os.path.join( fasta_path , "%s_%sJ.fasta"%(species, chain_type) ), read_all=True) 
            

    valignments = format_v_genes(valignments)
    jalignments = format_j_genes(jalignments)

    #ccalignments = format_c_genes(ccalignments, 'CC')
    #c1alignments = format_c_genes(c1alignments, 'C1')
    #c2alignments = format_c_genes(c2alignments, 'C2')
    #c3alignments = format_c_genes(c3alignments, 'C3')


    all_valignments = format_v_genes(all_valignments)
    all_jalignments = format_j_genes(all_jalignments)

    #all_C_alignments = { "CC":ccalignments,"C1":c1alignments,"C2":c2alignments,"C2":c2alignments}
    
    write_germlines( all_valignments, all_jalignments )

    # Combine the alignments to make putative germline alignments (obviously no d gene in there for Hs)
    # Write them to a stockholm alignment file.    
    combined_sequences  = make_putative_alignments( valignments, jalignments )

    # Write the constant domains each to file.
    #output_C_alignments(ccalignments, 'CC')
    #output_C_alignments(c1alignments, 'C1')
    #output_C_alignments(c2alignments, 'C2')
    #output_C_alignments(c3alignments, 'C3')


main()



