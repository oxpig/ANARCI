# Run the pipeline to get the sequences, format them and build the databases

DIR="."

# Rip the sequences from the imgt website. HTML may change in the future. 
mkdir -p $DIR/IMGT_sequence_files/htmlfiles
mkdir -p $DIR/IMGT_sequence_files/fastafiles
python3 $DIR/RipIMGT.py

# Format the alignments and handle imgt oddities to put into a consistent alignment format
mkdir -p $DIR/curated_alignments
mkdir -p $DIR/muscle_alignments
python3 $DIR/FormatAlignments.py

# Build the hmms for each species and chain.
# --hand option required otherwise it will delete columns that are mainly gaps. We want 128 columns otherwise ARNACI will fall over.
mkdir -p $DIR/HMMs
python3 $DIR/BuildHMM.py $DIR/HMMs/ALL.hmm $DIR/curated_alignments/ALL.stockholm
