import pytest
from anarci import run_pyhmmer


def test_run_pyhmmer():
    sequence = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS'
    hit_table, state_vectors, top_descriptions = run_pyhmmer(
        [('id', sequence)],
        hmm_database='ALL',
        ncpu=None,
        bit_score_threshold=80,
        hmmer_species=['human', 'mouse']
    )[0]

    state_vector = state_vectors[0]

    hmm_states = [hmm_state for (hmm_state, state_type), sequence_index in state_vector]
    state_types = ''.join(state_type for (hmm_state, state_type), sequence_index in state_vector)
    sequence_reproduced = ''.join(sequence[sequence_index] if sequence_index is not None else '' for _, sequence_index in state_vector)

    assert sequence == sequence_reproduced
    #                      QVQLQQSGA-ELARPGASVKMSCKASGYTF----TRYTMHWVKQRPGQGLEWIGYINPS--RGYTNYNQKFK-DKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS
    assert state_types == 'mmmmmmmmmdmmmmmmmmmmmmmmmmmmmmddddmmmmmmmmmmmmmmmmmmmmmmmmmddmmmmmmmmmmmdmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmiiiimmmmmmmmmmmmmmmmmm'
    assert hmm_states == list(range(1, 111)) + [111, 111, 111, 111, 111] + list(range(112, 129))

    assert len(hit_table) == 3
    assert hit_table[0] == ['id', 'description', 'evalue', 'bitscore', 'bias', 'query_start', 'query_end']
    assert hit_table[1][0] == 'mouse_H'
    assert hit_table[1][1] == ''
    assert hit_table[1][2] == pytest.approx(1.6e-52)
    assert hit_table[1][3] == pytest.approx(188.1, rel=1e-2)
    assert hit_table[1][4] == pytest.approx(3.3, rel=1e-2)
    assert hit_table[1][5] == 0
    assert hit_table[1][6] == 124

    assert len(top_descriptions) == 1
    assert top_descriptions[0]['id'] == 'mouse_H'
    assert top_descriptions[0]['description'] == ''
    assert top_descriptions[0]['evalue'] == pytest.approx(1.6e-52)
    assert top_descriptions[0]['bitscore'] == pytest.approx(188.1, rel=1e-2)
    assert top_descriptions[0]['bias'] == pytest.approx(3.3, rel=1e-2)
    assert top_descriptions[0]['query_start'] == 0
    assert top_descriptions[0]['query_end'] == 124
    assert top_descriptions[0]['chain_type'] == 'H'
    assert top_descriptions[0]['species'] == 'mouse'