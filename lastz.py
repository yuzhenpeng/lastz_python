#!/bin/python

#import os
#import tempfile

from subprocess import call


class LastZ:
    def __init__(self, output_dir=".", score_matrix=None, query="",
                 target="", output_format="lav", output_file=None):
        self.output_dir = output_dir
        self.score_matrix = score_matrix
        self.query = query
        self.target = target
        self.output_format = output_format
        self.output_file = output_file
        self.lastz()

    def lastz(self):
        """
        Parameter explanations:
            C=0 (--nochain):
                - Skip the chaining phase, not performed by default.
            O=600:
                - Gap open penalty
            E=150:
                - Gap extend penalty
            H=0 (--inner):
                - Perform additional alignment between gapped alignment blocks,
                  not performed by default.
            K=4500:
                - Score threshold for the x-drop extension method.
                - HSPs scoring lower than this are discarded.
                - By default, the entropy adjustment is used.
            L=3000:
                - Threshold for gapped extension, alignments scoring lower than
                  this are dicarded.
                - Gapped extension performed by default, and alignment ends are
                  trimmed to the locations giving the maximum score.
            M=254:
                - Dynamically mask the target sequence by excluding any
                  positions that appear in too many alignments from further
                  consideration for seeds.
            T=2:
                - Seed requires a 19-bp word with matches in 12 specific
                  positions (1110100110010101111).
            Y=15000:
                - Threshold for teminating gap extension. This restricts
                  the endpoints of each local alignment by limiting the local
                  region around each anchor in which extension is performed.
        """
        out_format = "--format=" + str(self.output_format)
        out_file = "--output=" + str(self.output_file)
        score_mat = "Q=" + str(self.score_matrix)

        call(["/nfs/research1/goldman/conor/tools/lastz/src/lastz", self.query,
              self.target, "C=0", "E=150", "H=0", "K=4500", "L=3000", "M=254",
              "O=600", "T=2", "Y=15000", score_mat, "--markend",
              out_format, out_file])

    #@staticmethod
    #def create_matrix_file(my_scoring_matrix):
    #    with tempfile.TemporaryFile() as tmp:
    #        tmp.write(my_scoring_matrix)


def main():
    with open("human_matrix.txt", "r") as f:
        scoring_matrix = f.read()
 
    print scoring_matrix

    LastZ(query="seq1.fa", target="seq2.fa", score_matrix="human_matrix.txt",
         output_file="out.fa")

    #LastZ.create_matrix_file(human_matrix)

    #LastZ.lastz()
    #
    #os.umask(saved_umask)


if __name__ == "__main__":
    main()
