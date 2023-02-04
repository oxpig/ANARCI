"""
Build All.HMM from stockholm aligned file using PyHMMER
"""

import sys
import pyhmmer


def build(out_path, in_path):
    alphabet = pyhmmer.easel.Alphabet.amino()
    print('Loading MSA from {}'.format(in_path))
    with pyhmmer.easel.MSAFile(in_path, digital=True, alphabet=alphabet) as msa_file:
        msas = list(msa_file)

    builder = pyhmmer.plan7.Builder(alphabet, architecture='hand')
    background = pyhmmer.plan7.Background(alphabet)

    print('Building HMMs from {} MSAs'.format(len(msas)))
    hmms = []
    with open(out_path, "wb") as output_file:
        for msa in msas:
            hmm, _, _ = builder.build_msa(msa, background)
            hmm.write(output_file)
            hmms.append(hmm)
    print('Pressing HMMs')
    pyhmmer.hmmpress(hmms, out_path)
    print('Saved HMMs to {}'.format(out_path))


if __name__ == '__main__':
    build(sys.argv[1], sys.argv[2])