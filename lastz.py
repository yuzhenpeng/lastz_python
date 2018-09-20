#!/bin/python

import datetime

from subprocess import call
from subprocess import check_output
from subprocess import STDOUT
from subprocess import PIPE
from subprocess import Popen
from sys import exit
from sys import argv
from os import remove
from os import devnull
from os import path
from os import makedirs


class LastZ:
    def __init__(self, output_dir=".", score_matrix=None, query="",
                 target="", output_format="lav", output_files=None,
                 log_path=None):
        """
        Python wrapper for performing a human-human lastz alignment.
        """
        # Declare variables
        self.temp_file_dir = \
            ("/hps/nobackup/research/goldmans/conor/lastz_temp_files/") + \
            query.split("/")[-1].split(".")[0] + "/"
        self.temp_file = self.temp_file_dir + \
            query.split("/")[-1].split(".")[0] + "." + \
            query.split("/")[-1].split(".")[1]
        self.FNULL = open(devnull, "w")
        self.log_file = log_path + "/" + str(target).split("/")[-1][:-5] + \
            "_" + str(query).split("/")[-1][:-5] + ".log"
        self.output_dir = output_dir
        self.score_matrix = score_matrix
        self.target = target
        self.target_only = target.split("/")[-1].split(".")[0]
        self.target_sizes = str(target)[:-5] + ".sizes"
        self.query = query
        self.query_only = query.split("/")[-1].split(".")[0]
        self.query_sizes = str(query)[:-5] + ".sizes"
        self.output_format = output_format
        self.output_file = self.temp_file + "_temp.lav"
        self.out_psl_file = self.temp_file + "_temp.psl"
        self.out_chain_file = self.temp_file + "_temp.chain"
        self.out_chainmergesort_file = self.temp_file + "_temp.chainMergeSort"
        self.all_chain_file = self.temp_file + ".allchain"
        self.pre_net_file = self.temp_file + ".prenet"
        self.target_net_file = self.temp_file + ".net"
        self.query_net_file = self.temp_file + ".net"
        self.syntenic_file = self.temp_file_dir + ".syntenicnet"
        self.unsorted_axt_file = self.temp_file_dir + "_unsorted.axt"
        self.axt_file = "aligned/" + self.query_only + "/axt/" + \
            self.target_only + "_" + self.query_only + ".axt"
        self.maf_file = "aligned/" + self.query_only + "/maf/" + \
            self.target_only + "_" + self.query_only + ".maf"

        # Call functions to proceed with alignment steps
        self.create_temp_sample_dir()
        self.check_2bit()
        self.sizes_files()
        self.lastz()
        self.validate_completion()
        self.convert_lav_to_psl()
        self.chaining()
        self.chain_merge_sort()
        self.chain_pre_net()
        self.chain_net()
        self.net_syntenic()
        self.net_to_axt()
        self.axt_sort()
        self.axt_to_maf()
        self.clean_up()

    def create_temp_sample_dir(self):
        if not path.exists(self.temp_file_dir):
            makedirs(self.temp_file_dir)

    def check_2bit(self):
        """
        Check the input sequence is in .2bit format.
        """
        print self.temp_file_dir
        print self.temp_file
        print self.output_file
        if self.target[-4:] != "2bit" or self.query[-4:] != "2bit":
            print "Error: input files must be in .2bit format."
            exit()

    def update_log(self, step_name, step_no, time):
        """
        Appends the alignment status to an associated log file.
        """
        # Remove milliseconds from the logged time
        print "log:", self.log_file
        time = str(time).split(".")[0]
        with open(self.log_file, "a+") as log_file:
            log_file.write("[{2}] {0}/12: {1}\n".format(
                step_no, step_name, time))
        print("[{2}] {0}/12: {1}".format(step_no, step_name, time))

    def sizes_files(self):
        """
        Check if the .sizes files exist, required later in the pipeline.
        """
        # Check files exist or not
        target_sizes = self.target[:-5] + ".sizes"
        query_sizes = self.query[:-5] + ".sizes"
        if not path.exists(target_sizes):
            p1 = Popen(["/homes/cwalker/tools/lastz/data/bin/twoBitInfo",
                       self.target, "stdout"], stdout=PIPE)
            call(["sort", "-k2rn"], stdin=p1.stdout,
                 stdout=open(target_sizes, "w"))
        if not path.exists(query_sizes):
            p1 = Popen(["/homes/cwalker/tools/lastz/data/bin/twoBitInfo",
                       self.query, "stdout"], stdout=PIPE)
            call(["sort", "-k2rn"], stdin=p1.stdout,
                 stdout=open(query_sizes, "w"))

    def lastz(self):
        """
        Parameter explanations:
            --nochain:
                - Skip the chaining phase, not performed by default.
            O=400:
                - Gap open penalty
            E=30:
                - Gap extend penalty
            H=3000 (--inner):
                - Perform additional alignment between gapped alignment blocks,
                  not performed by default.
            K=5000:
                - Score threshold for the x-drop extension method.
                - HSPs scoring lower than this are discarded.
                - By default, the entropy adjustment is used.
            L=5000:
                - Threshold for gapped extension, alignments scoring lower than
                  this are dicarded.
                - Gapped extension performed by default, and alignment ends are
                  trimmed to the locations giving the maximum score.
            M=10:
                - Dynamically mask the target sequence by excluding any
                  positions that appear in too many alignments from further
                  consideration for seeds.
            T=1:
                - Seed requires a 19-bp word with matches in 12 specific
                  positions (1110100110010101111).
            Y=15000:
                - Threshold for teminating gap extension. This restricts
                  the endpoints of each local alignment by limiting the local
                  region around each anchor in which extension is performed.
        """
        # Update the status log
        self.update_log("LASTZ alignment", "1", datetime.datetime.now())

        # Format options to pass to lastz
        out_format = "--format=" + str(self.output_format)
        out_file = "--output=" + str(self.output_file)
        score_mat = "Q=" + str(self.score_matrix)
        
         

        query_range = self.query + "[1..25000000]"
        target_range = self.target + "[1..25000000]"

        call(["/nfs/research1/goldman/conor/tools/lastz/src/lastz", target_range,
              query_range, "--nochain", "E=30", "H=3000", "K=5000", "L=5000",
              "M=10", "O=400", "T=1", "Y=15000", score_mat, "--markend",
              "--allocate:traceback=500.0M", out_format, out_file])

        print "After call"

    def validate_completion(self):
        self.update_log("Validating alignment", "2", datetime.datetime.now())
        """
        Check if lastz alignment completed successfully.
        """
        # Get last line of .lav file produced by LastZ.lastz()
        last_line = str(check_output(["tail", "-n", "1",
                                      self.output_file]).rstrip())
        # If last line is not the completion string, exit
        if last_line != "# lastz end-of-file":
            self.update_log("FAILED", "1", datetime.datetime.now())
            print("Error: {0} and {1} failed to align.".format(self.target,
                                                               self.query))
            exit()

    def convert_lav_to_psl(self):
        """
        Convert the generated lav file to psl format. Required for next axtChain.
        """
        self.update_log("Converting the .lav to .psl", "3",
                        datetime.datetime.now())
        # Name for psl file
        # Call kentUtils tool lavToPsl to perform the format conversion
        call(["/homes/cwalker/tools/lastz/data/bin/lavToPsl", self.output_file,
              self.out_psl_file], stdout=self.FNULL, stderr=STDOUT)
        # Delete the lav file
        remove(self.output_file)

    def chaining(self):
        """
        Wrapper for KentUtils axtChain.
        If two adjacent alignments are similar enough, they are concatenated
        into one fragment. These chain files are then sorted and combined into
        a single file.
        """
        self.update_log("Initial chaining", "4", datetime.datetime.now())
        # Reads in previously used scoring matrix
        score_scheme = "-scoreScheme=" + str(self.score_matrix)
        # Set to medium because closely related (same species...)
        linear_gap = "-linearGap=medium"
        # Set minScore to 5000
        min_score = "-minScore=5000"
        # Call axtChain
        call(["/homes/cwalker/tools/lastz/data/bin/axtChain", "-psl",
              min_score, score_scheme, linear_gap,
              self.out_psl_file, self.target, self.query,
              self.out_chain_file], stdout=self.FNULL, stderr=STDOUT)
        # Delete psl file
        remove(self.out_psl_file)

    def chain_merge_sort(self):
        """
        Perform chain, merge, sort step (combines sorted chained files into
        larger sorted files).
        """
        # Call to chainMergeSort from kentUtils
        self.update_log("Chaining, merging & sorting", "5",
                        datetime.datetime.now())
        call(["/homes/cwalker/tools/lastz/data/bin/chainMergeSort",
              self.out_chain_file], stdout=open(self.all_chain_file, "w"))
        # Remove file from previous step  containing shorter sorted chains
        remove(self.out_chain_file)

    def chain_pre_net(self):
        """
        Removes chains that don't have a chance of being netted.
        """
        self.update_log("Pre-netting processing", "6", datetime.datetime.now())
        # Call chainPreNet from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/chainPreNet",
              self.all_chain_file, self.target_sizes, self.query_sizes,
              self.pre_net_file], stdout=self.FNULL, stderr=STDOUT)
        # Remove all_chain file from previous step
        remove(self.all_chain_file)

    def chain_net(self):
        """
        Make alignment nets out of the chains.
        """
        # Call to chainNet from kentUtils
        self.update_log("Generating alignment nets", "7",
                        datetime.datetime.now())
        call(["/homes/cwalker/tools/lastz/data/bin/chainNet",
              self.pre_net_file, self.target_sizes, self.query_sizes,
              self.target_net_file, self.query_net_file],
             stdout=self.FNULL, stderr=STDOUT)

    def net_syntenic(self):
        """
        Add synteny information to the nets.
        """
        self.update_log("Adding synteny information", "8",
                        datetime.datetime.now())
        # Call to netSyntenic from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/netSyntenic",
              self.target_net_file, self.syntenic_file],
             stdout=self.FNULL, stderr=STDOUT)

    def net_to_axt(self):
        """
        Convert net and chains to .axt format.
        """
        self.update_log("Converting to .axt format", "9",
                        datetime.datetime.now())
        # Call to netToAxt from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/netToAxt",
              self.syntenic_file, self.pre_net_file, self.target,
              self.query, self.unsorted_axt_file], stdout=self.FNULL,
             stderr=STDOUT)

    def axt_sort(self):
        """
        Sort the generated axt file(s).
        """
        self.update_log("Sorting the .axt file(s)", "10",
                        datetime.datetime.now())
        # Call to axtSort from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/axtSort",
              self.unsorted_axt_file, self.axt_file],
             stdout=self.FNULL, stderr=STDOUT)

    def axt_to_maf(self):
        """
        Converts the generated .axt file to .maf format.
        """
        self.update_log("Converting .axt to .maf", "11", datetime.datetime.now())
        # Call to axtToMaf from kentUtils
        call(["/homes/cwalker/tools/lastz/data/bin/axtToMaf", self.axt_file,
              self.target_sizes, self.query_sizes, self.maf_file])

    def clean_up(self):
        """
        Removes any files which are not the final alignment.
        """
        self.update_log("Cleaning up temporary files", "12",
                        datetime.datetime.now())
        remove(self.syntenic_file)
        remove(self.unsorted_axt_file)
        remove(self.target_net_file)
        remove(self.query_net_file)
        remove(self.pre_net_file)
        # Update log to indicate completion
        with open(self.log_file, "a") as log_file:
            log_file.write("\n# alignment-finished\n")


def main():
    # Create directory to store logs of runs if it doesn't exist
    if not path.exists("lastz_logs"):
        makedirs("lastz_logs")

    # Define path for present sample, check if it exists, if not, create it
    sample_path = "lastz_logs" + "/" + str(argv[2]).split("/")[-1].split(".")[0]
    if not path.exists(sample_path):
        makedirs(sample_path)

    print "pr:",argv[1].split("/")[-1]

    # Name for log file of run
    log_file = str(argv[1]).split("/")[-1][:-5] + "_" + \
        str(argv[2]).split("/")[-1][:-5] + ".log"

    # Path to check for exitence of identical log file
    log_file_path = sample_path + "/" + log_file

    # Remove log file if it already existed
    if path.exists(log_file_path):
        remove(str(log_file_path))

    # Perform alignment, chaining and netting
    LastZ(target=argv[1],
          query=argv[2],
          score_matrix="/nfs/research1/goldman/conor/scripts/lastz/human_matrix.txt",
          output_files="out",
          log_path=sample_path)


if __name__ == "__main__":
    main()
