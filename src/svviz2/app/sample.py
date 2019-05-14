import codecs
import json
import logging
import numpy
import os, errno
import pickle
import pysam
import sys

from svviz2.io.readstatistics import ReadStatistics
from svviz2.utility.bam import bam_sort_index
from svviz2.utility.misc import str_to_bool

logger = logging.getLogger(__name__)


def _get_bam_headers(variant, allele):
    seqs = variant.seqs(allele)
    header = {"HD": {"VN": 1.3, "SO": "unsorted"}}
    sq = []
    for name in seqs:
        sq.append({"SN": name.replace("/", "__"), "LN": len(seqs[name])})
    header["SQ"] = sq
    return header


def get_sequencer_from_bam_header(bam):
    sequencer = None
    try:
        header = bam.header.get("RG", [])
    except AttributeError:  # pysam >= 0.14
        header = bam.header.to_dict().get("RG", [])

    for rg in header:
        if "PL" in rg:
            sequencer = rg["PL"].lower()

    return sequencer


class Sample(object):
    def __init__(self, name, bam_path, outdir, datahub, extra_args=None):
        self.name = name
        self.datahub = datahub
        self.outdir = outdir
        self.single_ended = False
        self.orientations = None
        self._search_distance = None
        dir_bams = os.path.join(self.outdir, "bams")
        try:
            os.mkdir(dir_bams)
        except:
            print("{} already created".format(dir_bams))
            pass

        symbolic_bam = os.path.join(dir_bams, "%s.bam" % self.name)
        symbolic_bam_bai = os.path.join(dir_bams, "%s.bam.bai" % self.name)
        try:
            os.symlink(bam_path, symbolic_bam)
            os.symlink("%s.bai" % bam_path, symbolic_bam_bai)
        except:
            pass

        self.bam_path = symbolic_bam
        self._bam = None

        ###SYMLINK
        # self.symbolic_bam_index = os.path.join(self.symbolic_bam, ".bai")
        self.outbams = {}
        self.outbam_paths = {}

        self.read_filter = None

        self._load(extra_args)

    def _load(self, extra_args):
        if not os.path.exists(self.bam_path):
            raise FileNotFoundError("Could not find bam file {}".format(self.bam_path))
        ##symbolic link for avoiding permission problems

        stats_file_path = self.bam_path + ".svviz_stats"
        print(stats_file_path)

        try:
            with open(stats_file_path) as inf:
                data = json.load(inf)

            self.single_ended = data["single_ended"]
            self.sequencer = data["sequencer"]
            self.max_base_quality = data["max_base_quality"]

            # this is a little hairy -- ideally we'd json serialize
            # the ReadStatistics object rather than pickle it
            self.read_statistics = pickle.loads(
                codecs.decode(data["read_statistics"].encode(), "base64")
            )
            logger.info("Loaded pre-computed read statistics for {}".format(self.name))

        except:
            self._load_from_bam(extra_args)
            read_stats_data = codecs.encode(
                pickle.dumps(self.read_statistics), "base64"
            ).decode()

            result = {
                "single_ended": self.single_ended,
                "sequencer": self.sequencer,
                "read_statistics": read_stats_data,
                "max_base_quality": self.max_base_quality,
            }

            with open(stats_file_path, "w") as outf:
                json.dump(result, outf)

    def _load_from_bam(self, extra_args):
        logger.info(
            "Calculating read statistics for {}".format(os.path.basename(self.bam_path))
        )

        import time

        t0 = time.time()
        self.read_statistics = ReadStatistics(self.bam)
        t1 = time.time()
        logger.info("TIME to get read statistics:{:.1f}s".format(t1 - t0))

        if self.read_statistics.orientations == "any":
            self.single_ended = True

        self.sequencer = "illumina"
        self.max_base_quality = 40.0
        if self.single_ended:
            mismatches = numpy.mean(self.read_statistics.number_mismatches)
            lengths = numpy.mean(self.read_statistics.readLengths)
            mismatch_rate = mismatches / lengths

            # these error rates aren't really accurate in terms of describing the
            # two platforms anymore, but they do correspond to the presets that
            # bwa mem has, which we're mimicking
            if mismatch_rate > 0.10:
                self.sequencer = "nanopore"
            elif mismatch_rate > 0.01:
                self.sequencer = "pacbio"
            elif numpy.isnan(mismatches) and lengths > 1000:
                self.sequencer = "pacbio"

        sequencer = get_sequencer_from_bam_header(self.bam)
        if sequencer in ["illumina", "pacbio", "nanopore"]:
            self.sequencer = sequencer

        if "sequencer" in extra_args:
            self.sequencer = extra_args["sequencer"].lower()
            assert self.sequencer in ["illumina", "pacbio", "nanopore"]
            self.single_ended = True
        if "single_ended" in extra_args:
            self.single_ended = str_to_bool(extra_args["single_ended"])
        if "max_base_quality" in extra_args:
            self.max_base_quality = float(extra_args["max_base_quality"])

        logger.info("Using realignment presets for {}".format(self.sequencer))

    @property
    def search_distance(self):
        if self._search_distance is None:
            if self.single_ended:
                self._search_distance = 1000
                if self.datahub.args.fast:
                    self._search_distance = 150
            else:
                search_distance = numpy.percentile(self.read_statistics.insertSizes, 99)
                self._search_distance = int(search_distance)
        return self._search_distance

    @property
    def align_distance(self):
        if self.single_ended:
            longest_reads = numpy.percentile(self.read_statistics.readLengths, 95) * 1.5
            return int(longest_reads)
        else:
            longest_inserts = (
                numpy.percentile(self.read_statistics.insertSizes, 98) * 1.5
            )
            return int(longest_inserts)

    @property
    def bam(self):
        if self._bam is None:
            self._bam = pysam.AlignmentFile(self.bam_path, "rb")
            try:
                self._bam.fetch()
            except ValueError:
                logger.error(
                    "ERROR: Need to create index for input bam file: {}".format(
                        self.bam_path
                    )
                )
                sys.exit(0)
        return self._bam

    def add_realignments(self, aln_sets):
        for allele in ["alt", "ref", "amb"]:
            self.outbam(allele, "w")

        for aln_set in aln_sets:
            # if aln_set.supports_allele != "amb":
            if aln_set.supporting_aln is None:
                continue
            aln_set.supporting_aln.fix_flags()

            outbam = self.outbam(aln_set.supports_allele, "w")
            if self.single_ended:
                if aln_set.supports_allele == "amb":
                    if aln_set.supporting_aln.chrom not in self.datahub.variant.seqs(
                        "amb"
                    ):
                        continue

                outbam.write(aln_set.supporting_aln._read)
            else:
                if aln_set.supports_allele == "amb":
                    if (
                        aln_set.supporting_aln.aln1.chrom
                        not in self.datahub.variant.seqs("amb")
                    ):
                        continue
                    if (
                        aln_set.supporting_aln.aln2.chrom
                        not in self.datahub.variant.seqs("amb")
                    ):
                        continue

                outbam.write(aln_set.supporting_aln.aln1._read)
                outbam.write(aln_set.supporting_aln.aln2._read)

    def finish_writing_realignments(self):
        for allele in ["alt", "ref", "amb"]:
            self.outbams[allele].close()
            self.outbams.pop(allele)
            try:
                bam_sort_index(self.outbam_paths[allele])
            except:
                print("ERROR!" * 30)
                raise

    def has_realignments(self):
        for allele in ["alt", "ref", "amb"]:
            try:
                self.outbam(allele, "r")
            except OSError:
                return False
        return True

    def outbam(self, allele, mode):
        if mode == "w":
            if not allele in self.outbams:
                self.outbams[allele] = pysam.AlignmentFile(
                    self.outbam_paths[allele],
                    "wb",
                    header=_get_bam_headers(self.datahub.variant, allele),
                )
            return self.outbams[allele]

        assert not allele in self.outbams, "forgot to close outbam before re-opening"

        return pysam.AlignmentFile(
            self.outbam_paths[allele].replace(".bam", ".sorted.bam")
        )

    def __getstate__(self):
        """ allows pickling of Samples()s """
        state = self.__dict__.copy()
        del state["bam"]
        return state

