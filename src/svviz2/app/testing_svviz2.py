import logging
import sys
import time
from argparse import Namespace
import collections
import json
import numpy
import os
import pandas


from svviz2.app import commandline
from svviz2.app.datahub import DataHub
from svviz2.visualize import visualize

# from svviz2.app import report
from svviz2.visualize import dotplots

# report
from svviz2.app.variants import non_negative
from svviz2.utility import statistics
from svviz2.remap import genotyping


FORMAT = "%(asctime)s - %(name)-25s - %(levelname)-5s - %(message)s"
DATEFMT = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
logger = logging.getLogger(__name__)


def report(datahub):
    results = []
    results.extend(tally_support(datahub))
    # results.extend(tally_segments(datahub))
    # results.extend(tally_nearby_polymorphisms(datahub))

    result_df = pandas.DataFrame(results, columns=["sample", "allele", "key", "value"])
    result_df["event"] = datahub.variant.name
    print(result_df)
    report_filename = "{}.report.tsv".format(datahub.variant.short_name())
    report_path = os.path.join(datahub.args.outdir, report_filename)
    result_df.to_csv(report_path, sep="\t", index=False)
    return report_filename


def tally_support(datahub):
    results = []
    for sample_name, sample in datahub.samples.items():
        allele_count = {}
        allele_mapq_sum = {}
        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")
            cur_results, count, mapq_sum, _ = _tally_support(bam)
            allele_count[allele] = count
            allele_mapq_sum[allele] = mapq_sum
            for cur_result in cur_results:
                results.append((sample_name, allele) + cur_result)
    return results


def _tally_support(bam):
    count = 0
    weighted_count = 0
    mapq_sum = 0
    breakpoint_overlaps = collections.defaultdict(list)
    breakpoint_counts = collections.Counter()
    extensions = collections.defaultdict(list)

    breakpoints = set()
    for read in bam:
        if not read.is_paired or read.is_read1:
            count += 1
            weighted_count += read.mapq / 40.0
            mapq_sum += 1 - statistics.phred_to_prob(read.mapq, 10)
            cur_breakpoint_overlaps = json.loads(read.get_tag("Ov"))
            for breakpoint, info in cur_breakpoint_overlaps.items():
                breakpoints.add(breakpoint)
                overlap, overlaps_sequence, extension = info
                breakpoint_overlaps[breakpoint].append(overlap)
                breakpoint_counts[(breakpoint, overlaps_sequence)] += 1
                extensions[breakpoint].append(extension)
    results = [("count", count)]
    return results, count, mapq_sum, weighted_count


def get_datahub():
    args = Namespace(
        align_distance=None,
        aligner="bwa",
        also_plot_context=None,
        bam=[
            "/home/varelad/data/I-H-134709-N1-1-D1-1.bam",
            "/home/varelad/data/I-H-134709-T1-1-D1-1.bam",
            # "/home/varelad/data/I-H-134709-T2-1-D1-1.bam",
            # "/home/varelad/data/I-H-134709-T3-1-D1-1.bam",
            # "/home/varelad/data/I-H-134709-T4-1-D1-1.bam",
        ],
        batch_size=10000,
        dotplots_only=False,
        downsample=None,
        fast=False,
        first_variant=None,
        format="pdf",
        last_variant=None,
        min_mapq=None,
        no_dotplots=True,
        no_render=True,
        no_report=False,
        only_plot_context=None,
        only_realign_locally=False,
        outdir="/home/varelad/testing",
        ref="/work/isabl/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta",
        render_only=False,
        report_only=True,
        savereads=False,
        variants="/home/varelad/subset_709.vcf",
    )
    datahub = DataHub()
    datahub.set_args(args)
    datahub.align_distance = 0

    for _, sample in datahub.samples.items():
        logger.info("Search distance: {:,}bp".format(sample.search_distance))

    datahub.align_distance = max(
        sample.align_distance for sample in datahub.samples.values()
    )
    if datahub.args.align_distance is not None:
        assert (
            datahub.args.align_distance > 0
        ), "--align-distance must be a positive integer"
        datahub.align_distance = datahub.args.align_distance
    logger.info("Align distance: {:,}bp".format(sample.align_distance))
    return datahub


def run(datahub):
    """ this runs the app on the provided datahub """
    list_reports = []
    file_name_reports = os.path.join(datahub.args.outdir, "names_reports.txt")
    for _ in datahub.get_variants():
        if datahub.should_genotype:
            t0 = time.time()
            datahub.genotype_cur_variant()
            t1 = time.time()
            print("TIME:::", t1 - t0)
            # print(datahub.variant)
        if datahub.should_render:
            visualize.visualize(datahub)
        if datahub.should_generate_reports:
            name_file = report(datahub)
        list_reports.append(name_file)
        with open(file_name_reports, "w") as f:
            for item in list_reports:
                f.write("%s\n" % item)
        if datahub.should_generate_dotplots:
            dotplots.generate_dotplots(datahub)
    datahub.cleanup()


datahub = get_datahub()
run(datahub)
print("DONE")

