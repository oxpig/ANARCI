
# Import the anarci function.
from anarci import anarci

# Format the sequences that we want to number. 
sequences = [ ("12e8:H","EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAP"),
              ("12e8:L","DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASV"),
              ("scfv:A","DIQMTQSPSSLSASVGDRVTITCRTSGNIHNYLTWYQQKPGKAPQLLIYNAKTLADGVPSRFSGSGSGTQFTLTISSLQPEDFANYYCQHFWSLPFTFGQGTKVEIKRTGGGGSGGGGSGGGGSGGGGSEVQLVESGGGLVQPGGSLRLSCAASGFDFSRYDMSWVRQAPGKRLEWVAYISSGGGSTYFPDTVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARQNKKLTWFDYWGQGTLVTVSSHHHHHH"),
              ("lysozyme:A","KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL")]

# Hand the list of sequences to the anarci function. Number them with the IMGT scheme
results = anarci(sequences, scheme="imgt", output=False)

# Unpack the results. We get three lists
numbering, alignment_details, hit_tables = results

# Each has the same number of elements as the number of sequences submitted
assert len(numbering) == len(alignment_details) == len(hit_tables) == len( sequences )

print('I am using the anarci function to number and get all the details about the following sequences')
print(sequences)

print('\n')
# Iterate over the sequences
for i in range(len(sequences)):    
    if numbering[i] is None:
        print('ANARCI did not number', sequences[i][0])
    else:
        print('ANARCI numbered', sequences[i][0])
        print('It identified %d domain(s)'%len(numbering[i]))
        
        # Iterate over the domains
        for j in range(len(numbering[i])):
            domain_numbering, start_index, end_index = numbering[i][j]
            print('This is the IMGT numbering for the %d\'th domain:'%j, domain_numbering)
            print('This is the bit of the sequence it corresponds to:', sequences[i][1][start_index:end_index+1])
            print('These are the details of the alignment:')
            for (key,value) in alignment_details[i][j].items():
                print(key, ':', value)
            print('This is the summary of the hits that HMMER found')
            for line in hit_tables[i]:
                print(line)      
    print('\n','_'*40)
    
print('Do with this infomation as you wish')

print('\n','*'*40)

# Want to just get a quick numbering without caring about the other details?
from anarci import number

seq = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAP"
numbering, chain_type = number( seq, scheme = 'kabat' )
print('Alternatively we can simply number the first domain of a sequence with the number function')
print('I gave it this sequence\n', seq)
print('ANARCI told me it was a', chain_type, 'chain')
print('This is the first domain\'s Kabat numbering:')
print(numbering)
