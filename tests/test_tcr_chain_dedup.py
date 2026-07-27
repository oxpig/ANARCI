"""
Regression test for a domain-dedup bug that affects any chain type search
restriction (allow={'H'}, {'K','L'}, {'A','B'}, etc.), not just this example.

Before the fix, _domains_are_same() in anarci.py considered any two
overlapping hits to be the same domain regardless of chain type, so
whenever a weaker hit of an allowed type overlapped a stronger hit of a
disallowed type, the weaker hit was discarded before "allow" filtering
ever saw it, and the sequence was incorrectly left unnumbered.

This is easiest to reproduce with TCR alpha/delta (TRA/TRD), since their
V-regions are nearly identical and always overlap: a human TCR alpha chain
can score highest against the delta HMM, so restricting to allow={'A','B'}
incorrectly returned None even though a strong alpha hit was in the hit
table.
"""

from anarci import anarci

# TCR sequence whose top HMM hit is human_D (delta) but whose second-best hit
# is human_A (alpha). Restricting to allow={'A','B'} should recover the alpha hit
# rather than discarding it in the dedup step.
sequences = [("seq1", "AQKITQTQPGMFVQEKEAVTLDCTYDTSDPSYGLFWYKQPSSGEMIFLIYQGSYDQQNATEGRYSLNFQKARKSANLVISASQLGDSAMYFCAMSQLNSGNTPLVFGKGTRLSVIANIONPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPE")]

numbering, alignment_details, hit_tables = anarci(sequences, scheme="imgt", allow={"A", "B"}, output=False)

assert numbering[0] is not None, "Expected the alpha chain to be numbered when allow={'A','B'}, got None"
assert alignment_details[0][0]["chain_type"] == "A", "Expected the recovered domain to be chain type A (alpha)"

print("PASS: TRA sequence numbered correctly despite a stronger overlapping TRD hit")
