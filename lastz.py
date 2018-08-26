#!/bin/python

#import os
#import tempfile

from subprocess import call
from subprocess import check_output
from sys import exit
from os import remove


class LastZ:
    def __init__(self, output_dir=".", score_matrix=None, query="",
                 target="", output_format="lav", output_files=None):
        """
        Python wrapper for performing a human-human lastz alignment.
        """
        self.output_dir = output_dir
        self.score_matrix = score_matrix
        self.target = target
        self.target_sizes = str(target)[:-5] + ".sizes"
        self.query = query
        self.query_sizes = str(query)[:-5] + ".sizes"
        self.output_format = output_format
        self.output_file = output_files + ".lav"
        self.out_psl_file = output_files + ".psl"
        self.out_chain_file = output_files + ".chain"
        self.out_chainmergesort_file = output_files + ".chainMergeSort"
        self.all_chain_file = str(target)[:-5] + "_" + \
            str(query)[:-5] + ".allchain"
        self.pre_net_file = str(target)[:-5] + "_" + str(query)[:-5] + ".prenet"
        self.target_net_file = str(target)[:-5] + ".net"
        self.query_net_file = str(query)[:-5] + ".net"
        self.syntenic_file = str(target)[:-5] + "_" + str(query)[:-5] + \
            ".syntenicnet"
        self.lastz()
        self.validate_completion()
        self.convert_lav_to_psl()
        self.chaining()
        self.chain_merge_sort()
        self.chain_pre_net()
        self.chain_net()
        self.net_syntenic()

    def lastz(self):
        """
        Parameter explanations:
            --nochain:
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

        call(["/nfs/research1/goldman/conor/tools/lastz/src/lastz", self.target,
              self.query, "--nochain", "E=150", "H=0", "K=4500", "L=3000",
              "M=254", "O=600", "T=2", "Y=15000", score_mat, "--markend",
              out_format, out_file])

    def validate_completion(self):
        """
        Check if lastz alignment completed successfully.
        """
        # Get last line of .lav file produced by LastZ.lastz()
        last_line = str(check_output(["tail", "-n", "1",
                                      self.output_file]).rstrip())
        # If last line is not the completion string, exit
        if last_line != "# lastz end-of-file":
            print("Error: {0} and {1} failed to align.".format(self.target,
                                                               self.query))
            exit()

    def convert_lav_to_psl(self):
        """
        Conver the generated lav file to psl format. Required for next axtChain.
        """
        # Name for psl file
        # Call kentUtils tool lavToPsl to perform the format conversion
        call(["/homes/cwalker/tools/lastz/data/bin/lavToPsl", self.output_file,
              self.out_psl_file])
        # Delete the lav file
        remove(self.output_file)

    def chaining(self):
        """
        Wrapper for KentUtils axtChain.
        If two adjacent alignments are similar enough, they are concatenated
        into one fragment. These chain files are then sorted and combined into
        a single file.
        """
        # Reads in previously used scoring matrix
        score_scheme = "-scoreScheme=" + str(self.score_matrix)
        # Set to medium because closely related (same species...)
        linear_gap = "-linearGap=medium"
        # Set minScore to 5000
        min_score = "-minScore=5000"
        # Call axtChain
        call(["/homes/cwalker/tools/lastz/data/bin/axtChain", "-psl",
              min_score, score_scheme, linear_gap,
              self.out_psl_file, self.target, self.query, self.out_chain_file])
        # Delete psl file
        remove(self.out_psl_file)

    def chain_merge_sort(self):
        """
        Perform chain, merge sort step (combines sorted chained files into
        larger sorted files).
        """
        # Call to chainMergeSort from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/chainMergeSort",
              self.out_chain_file], stdout=open(self.all_chain_file, "w"))
        # Remove file from previous step  containing shorter sorted chains
        remove(self.out_chain_file)

    def chain_pre_net(self):
        """
        Removes chains that don't have a chance of being netted.
        """
        # Call chainPreNet from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/chainPreNet",
              self.all_chain_file, self.target_sizes, self.query_sizes,
              self.pre_net_file])
        # Remove all_chain file from previous step
        remove(self.all_chain_file)

    def chain_net(self):
        """
        Make alignment nets out of the chains.
        """
        # Call to chainNet from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/chainNet",
              self.pre_net_file, self.target_sizes, self.query_sizes,
              self.target_net_file, self.query_net_file])

    def net_syntenic(self):
        """
        Add synteny information to the nets.
        """
        # Call to netSyntenic from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/netSyntenic",
              self.target_net_file, self.syntenic_file])


def main():
    LastZ(target="seq1.2bit", query="seq2.2bit",
          score_matrix="human_matrix.txt", output_files="out")


if __name__ == "__main__":
    main()
