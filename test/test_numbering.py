from anarci import number


def test_imgt_heavy():
    sequence = "QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS"
    numbering, chain_type = number(sequence, "imgt")

    assert chain_type == 'H'

    positions = [f'{pos}{letter}'.strip() for (pos, letter), aa in numbering]
    assert positions == [str(i) for i in range(1, 111)] + ['111', '111A', '111B', '112B', '112A'] + [str(i) for i in range(112, 129)]

    aligned_sequence = ''.join(aa for (pos, letter), aa in numbering)
    assert aligned_sequence == 'QVQLQQSGA-ELARPGASVKMSCKASGYTF----TRYTMHWVKQRPGQGLEWIGYINPS--RGYTNYNQKFK-DKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVSS'

    reproduced_sequence = ''.join(aa for (pos, letter), aa in numbering if aa != '-')
    assert reproduced_sequence == sequence


def test_imgt_light():
    sequence = "DVVMTQSPLSLPVTLGQPASISCRSSQSLVYSDGNTYLNWFQQRPGQSPRRLIYKVSNRDSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCMQGTFGQGTKVEIK"
    numbering, chain_type = number(sequence, "imgt")

    assert chain_type == 'L'  # L for light

    positions = [f'{pos}{letter}'.strip() for (pos, letter), aa in numbering]
    assert positions == [str(i) for i in range(1, 128)]

    aligned_sequence = ''.join(aa for (pos, letter), aa in numbering)
    assert aligned_sequence == 'DVVMTQSPLSLPVTLGQPASISCRSSQSLVYS-DGNTYLNWFQQRPGQSPRRLIYKV-------SNRDSGVP-DRFSGSG--SGTDFTLKISRVEAEDVGVYYCMQ---------GTFGQGTKVEIK'

    reproduced_sequence = ''.join(aa for (pos, letter), aa in numbering if aa != '-')
    assert reproduced_sequence == sequence


def test_imgt_heavy_incomplete():
    """Check that a few mismatching residues at N and C termini are treated as part of the aligned domain"""
    sequence = "AAQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVAA"
    numbering, chain_type = number(sequence, "imgt")

    aligned_sequence = ''.join(aa for (pos, letter), aa in numbering)
    assert aligned_sequence == 'AAQLQQSGA-ELARPGASVKMSCKASGYTF----TRYTMHWVKQRPGQGLEWIGYINPS--RGYTNYNQKFK-DKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYELVISCLDYWGQGTTLTVAA'

