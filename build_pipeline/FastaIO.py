#!/usr/bin/env python

"""
Custom class to parse FASTA files.
@author J. Leem
jinwoo.leem [~at~] gmail.com
"""

def chunkify(fh):
    """ Divides FASTA files into chunks. """
    chunk = ""
    chunk_start = False
    description, seq = "", ""

    for line in fh:
        line = line.strip()

        # Skip over empty lines.
        if not line:
            continue

        # FASTA records always start with ">"
        # Case 1: No chunk has ever been defined
        #         so we start a new one.
        if line[0] == ">" and not chunk_start:
            chunk_start = True
            description = line[1:]
            continue

        # Case 2: We have found at least one chunk,
        #         so let's yield a FASTARecord object when we see the next ">"
        elif line[0] == ">" and chunk_start:
            yield FASTARecord( description, chunk )

            # Now it's a new record, so start a new description
            # Wipe out the existing sequence.
            description = line[1:]
            chunk = ""

        elif line[0] != ">" and chunk_start:
            chunk += line

    yield FASTARecord( description, chunk )

class FASTARecord(object):
    def __init__(self, description, seq):
        self.description = description.strip(">")
        self.seq = seq
