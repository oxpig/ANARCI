```
ANARCI                                                 \\    //
Antibody Numbering and Antigen Receptor ClassIfication  \\  //
                                                          ||
(c) Oxford Protein Informatics Group (OPIG). 2015-20      ||

Author: James Dunbar (dunbar@stats.ox.ac.uk)
        Charlotte Deane (deane@stats.ox.ac.uk)

Contact: opig@stats.ox.ac.uk

```

# Usage:

* Numbering a single sequence
```python
ANARCI -i EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA
```

* Numbering sequences in a FASTA file

```python
ANARCI -i myfile.fasta 
```

* Please note that while ANARCI uses alignment to species V and J germlines to determine the species of the antibody for purposes of numbering, we do not recommend using ANARCI as your primary species annotation tool

# Installation

The easiest way to install ANARCI and its dependencies is using conda

```python
conda install -c conda-forge biopython -y
conda install -c bioconda hmmer=3.3.2 -y
cd ANARCI
python setup.py install
```

# Further info

## Output files

* The numbering file.

The numbering file (`--outfile` or stdout) reports the numbering for all sequences given in the sequence file. Each record is separated by "//".
Those chains for which no significant alignment was found report the name as in the fasta file. e.g:

```    
# 1A14:N|PDBID|CHAIN|SEQUENCE
//
     
    Those sequences where a significant alignment has been found report as below: 
    
# 1A14:H|PDBID|CHAIN|SEQUENCE
# ANARCI numbered
# Domain 1 of 1
# Most significant HMM hit
#|species|chain_type|e-value|score|seqstart_index|seqend_index|
#|mouse|H|8.6e-58|184.9|0|119|
# Scheme = imgt
H 1       Q
H 2       V
H 3       Q
H 4       L
H 5       Q
  . 
  .
  .
//

    where:
    species          = The species of the most significant aligned HMM
    chain_type       = The chain type of the most significant aligned HMM
    e-value          = The e-value of the alignment to the most significant aligned HMM
    score            = The bit-score of the alignment to the most significant aligned HMM
    seqstart_index   = The index in the sequence of the first numbered residue
    seqend_index     = The index in the sequence of the last numbered residue
    Scheme           = The numbering scheme used to number the sequence
    
    Then follows the numbering. Chain type (H, L (for both kappa(K) and lambda(L) chain types) , A (alpha), B (beta))

    If the "assign_germline" option has been specified the further following lines are added to the header. e.g.

# Most sequence-identical germlines
#|species|v_gene|v_identity|j_gene|j_identity|
#|mouse|IGHV1-12*01|0.86|IGHJ2*01|0.79|

    where:
    species          = The species of the most sequence identical germline
    v_gene           = The identifier of the most sequence identical germline over the v-region
    v_identity       = The sequence identity over the v-region to the most sequence identical germline
    j_gene           = The identifier of the most sequence identical germline over the j-region
    j_identity       = The sequence identity over the j-region to the most sequence identical germline
```

* The csv format output file.

    When the `--csv` option is specified, numbered sequences are output into separate comma separated value files depending on their
    chain type. This provides a horizontal output format and contains all the properties detailed above. In addition, sequences 
    are aligned according to the numbering scheme. 


* The hit file.

    The hit file reports the statistics for all alignments to each HMM in the database even if the sequence was not numbered.
    Each record is separated by "//". 
    
    The corresponding hit table for the numbered entry above looks like:

```    
    """
        NAME     1a14_H mol:protein length:120  NC10 FV (HEAVY CHAIN)
    SEQUENCE QVQLQQSGAELVKPGASVRMSCKASGYTFTNYNMYWVKQSPGQGLEWIGIFYPGNGDTSYNQKFKDKATLT
    SEQUENCE ADKSSNTAYMQLSSLTSEDSAVYYCARSGGSYRYDGGFDYWGQGTTVTV
                   id       description            evalue          bitscore              bias   best_dom_evalue best_dom_bitscore     best_dom_bias    domain_exp_num    domain_obs_num
              mouse_H                             1.1e-57             184.5               1.5           1.3e-57             184.4               1.5               1.0                 1
              human_H                             7.8e-53             169.0               1.9           8.7e-53             168.8               1.9               1.0                 1
                rat_H                             4.7e-47             150.2               2.2           5.2e-47             150.0               2.2               1.0                 1
             rabbit_H                             3.7e-37             118.2               0.7             4e-37             118.1               0.7               1.0                 1
                pig_H                             1.5e-35             113.3               2.7           1.6e-35             113.1               2.7               1.0                 1
             rhesus_H                             4.4e-32             101.5               1.8           4.9e-32             101.4               1.8               1.0                 1
              mouse_B                             2.4e-19              60.6               0.7           2.6e-19              60.5               0.7               1.0                 1
              human_B                             4.2e-19              59.7               0.9           4.6e-19              59.5               0.9               1.0                 1
              mouse_A                             8.7e-19              58.5               1.1           9.6e-19              58.4               1.1               1.0                 1
              human_A                             1.7e-18              57.6               0.9           1.9e-18              57.5               0.9               1.0                 1
              mouse_D                             5.1e-17              53.3               0.7           5.9e-17              53.1               0.7               1.1                 1
             rhesus_L                             1.6e-16              51.7               2.7           1.9e-16              51.4               2.7               1.1                 1
              human_L                             1.7e-15              48.3               3.5             2e-15              48.0               3.5               1.1                 1
              human_D                             6.7e-15              46.1               0.2           7.4e-15              45.9               0.2               1.0                 1
             rhesus_K                             3.9e-13              40.6               1.7           5.1e-13              40.2               1.7               1.2                 1
              mouse_G                             4.1e-13              40.3               0.0           4.3e-13              40.2               0.0               1.0                 1
             rabbit_L                             6.1e-13              40.0               2.8           8.1e-13              39.6               2.8               1.2                 1
                rat_K                             3.9e-12              37.4               1.4           4.4e-12              37.2               1.4               1.1                 1
                pig_L                             4.2e-12              37.5               1.0           4.7e-12              37.3               1.0               1.1                 1
              mouse_K                             1.2e-11              35.7               2.6           1.3e-11              35.6               2.6               1.1                 1
              human_K                             2.2e-11              34.8               2.9           3.5e-11              34.2               2.9               1.4                 1
              mouse_L                             1.9e-10              31.8               2.2           3.4e-10              30.9               2.2               1.4                 1
                rat_L                             2.5e-10              31.7               1.2           2.9e-10              31.5               1.2               1.1                 1
                pig_K                             3.2e-10              31.1               1.9           4.5e-10              30.6               1.9               1.3                 1
              human_G                             2.9e-09              27.8               0.8           4.9e-09              27.1               0.8               1.4                 1
             rabbit_K                             2.5e-06              18.4               5.8           4.2e-06              17.7               5.8               1.4                 1
    //
    """

```

We therefore get a ranking of the alignments to each chain type. 

## Schemes:

    
* Currently implemented schemes: 
    IMGT
    Chothia (IGs only)
    Kabat (IGs only)
    Martin / Enhanced Chothia (IGs only)
    AHo 
    Wolfguy (IGs only)
        
* Currently recognisable species (chains):
    Human (heavy, kappa, lambda, alpha, beta)
    Mouse (heavy, kappa, lambda, alpha, beta)
    Rat (heavy, kappa, lambda)
    Rabbit (heavy, kappa, lambda)
    Pig (heavy, kappa, lambda)
    Rhesus Monkey (heavy, kappa)
    
Other species may still be numbered correctly and the chain type recognised but the species be incorrect. e.g. llama VHH.


* IMGT     - has 128 possible positions for *all* antigen receptor types. These are supposed to be structurally equivalent.
            In theory these are supposed to account for all possible positions. However, insertions are possible especially
            at CDR3. ANARCI gives the insertion codes as letters. Insertions at CDR3 are placed symmetrically about imgt
            positions 111 and 112. e.g. 111-ABCD DCBA-112. 
            
* Kabat    - is defined for heavy and light chain antibody chains only (in ANARCI). Positions in the two chain types are not
            equivalent. Insertions occur at specific positions and can occur in both the framework and the CDRs. They are
            annotated from A->Z. e.g 100ABCDEFGH 101.       

* Chothia  - is defined for heavy and light chain antibody chains only (in ANARCI). Numbering in the two chain types are not
            equivalent. Insertions occur at specific positions and can occur in both the framework and the CDRs. They are
            annotated from A->Z. e.g 100ABCDEFGH 101. It differs to the Kabat scheme by the position it places the insertions
            at CDRH1.   

* Martin   - is defined for heavy and light chain antibody chains only. Numbering in the two chain types are not equivalent. 
            Insertions occur at specific positions and can occur in both the framework and the CDRs. They are annotated from
            A->Z. e.g 100ABCDEFGH 101. It differs to the Chothia scheme by the position it places the certain insertions in
            the framework. It is also referred to as the enhanced Chothia scheme.

* AHo      - has 149 possible for *all* antigen receptor types. These are supposed to be structurally equivalent. The AHo
            scheme's large number of positions is supposed to account for all possible positions without the need for 
            specifying insertion positions. In ANARCI, extra residues in the framework may be represented by insertions
            although these are unlikely to occur in natural sequences.

* Wolfguy  - is defined for heavy and light antibody chains. Numbering in the two chain types are not equivalent. Different 
            regions of the domain are denoted by a range of numbers. Many possible positions in the CDRs mean that insertion
            codes are not required. In ANARCI, extra residues in the framework may be represented by insertions
            although these are unlikely to occur in natural sequences. The CDRs are numbered in an 'up' and 'down' direction.
            The annotations of CDRL1 is defined according to the canonical structure. In ANARCI this is recognised by taking
            a sequence similarity to hard coded sequence motifs for different lengths.

