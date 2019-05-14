import logging
import sys
import time
from argparse import Namespace
import collections
import json
import numpy
import os
import pandas as pd
import itertools as it
from pysam import VariantFile


from svviz2.app import commandline
from svviz2.app.datahub import DataHub

FORMAT = "%(asctime)s - %(name)-25s - %(levelname)-5s - %(message)s"
DATEFMT = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
logger = logging.getLogger(__name__)


def report(datahub):
    results = []
    results.extend(tally_support(datahub))
    print("RESULTS")
    print(results)
    #######
    args = [iter(results)] * 2
    results_grouped = list(it.zip_longest(*args, fillvalue=None))
    print(results_grouped)
    list_dict_samples = []
    for sample in results_grouped:
        if sample[0][0] == sample[1][0]:
            sample_name = sample[0][0]
            alt = sample[0][3]
            ref = sample[1][3]
            if int(alt + ref) != 0:
                vaf = round(int(alt) / int(alt + ref), 3)
            else:
                vaf = 0
            dict_results = {"sample": sample_name, "alt,ref": [alt, ref], "vaf": vaf}
        else:
            print("Error in results")
        list_dict_samples.append(dict_results)

    # result_df = pandas.DataFrame(results, columns=["sample", "allele", "key", "value"])
    result_df = pd.DataFrame(list_dict_samples, columns=["sample", "alt,ref", "vaf"])
    result_df["event"] = datahub.variant.name
    print(result_df)
    report_filename = "{}_report.tsv".format(datahub.variant.short_name())
    report_path = os.path.join(datahub.args.outdir, report_filename)
    result_df.to_csv(report_path, sep="\t", index=False)
    return report_path


def report_option2(datahub):
    results = []
    results.extend(tally_support(datahub))
    print(results)
    if results[0][0] == results[1][0]:
        sample_name = results[0][0]
    alt = results[0][3]
    ref = results[1][3]
    if int(alt + ref) != 0:
        vaf = round(int(alt) / int(alt + ref), 3)
    else:
        vaf = 0
    dict_results = {"sample": sample_name, "alt,ref": [alt, ref], "vaf": vaf}
    print(dict_results)
    return dict_results, sample_name


def tally_support(datahub):
    results = []
    for sample_name, sample in datahub.samples.items():
        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")
            cur_results = _tally_support(bam)
            for cur_result in cur_results:
                results.append((sample_name, allele) + cur_result)
    return results


def _tally_support(bam):
    count = 0
    for read in bam:
        if not read.is_paired or read.is_read1:
            count += 1
    results = [("count", count)]
    return results


def get_datahub():
    args = Namespace(
        align_distance=None,
        aligner="bwa",
        bam=[
            "/ifs/res/leukgen/local/opt/leukdc/data/workflows/25/37/42537/data/bam/I-H-134709-T2-1-D1-1.bam",
            # "/home/varelad/data/I-H-134709-T1-1-D1-1.bam",
            # "/home/varelad/data/I-H-134709-T2-1-D1-1.bam",
            # "/home/varelad/data/I-H-134709-T3-1-D1-1.bam",
            # "/home/varelad/data/I-H-134709-T4-1-D1-1.bam",
        ],
        batch_size=10000,
        downsample=None,
        fast=False,
        first_variant=None,
        last_variant=None,
        min_mapq=None,
        no_report=False,
        outdir="/home/varelad/testing",
        ref="/work/isabl/ref/homo_sapiens/GRCh37d5/genome/gr37.fasta",
        report=True,
        savereads=False,
        variants="/home/varelad/subset_709_copy.vcf",
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


###paralellization by variant
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
        if datahub.should_generate_reports:
            path = report(datahub)
        list_reports.append(path)
        with open(file_name_reports, "w") as f:
            for item in list_reports:
                f.write("%s\n" % item)
    datahub.cleanup()


###paralellization by sample
def run_option2(datahub):
    """ this runs the app on the provided datahub """
    list_reports = []
    for variant in datahub.get_variants():
        if datahub.should_genotype:
            t0 = time.time()
            datahub.genotype_cur_variant()
            t1 = time.time()
            print("TIME:::", t1 - t0)
        if datahub.should_generate_reports:
            result_variant, sample_name = report_option2(datahub)
            result_variant["event"] = variant.name
        list_reports.append(result_variant)
        print(list_reports)
    result_df = pd.DataFrame(
        list_reports, columns=["sample", "alt,ref", "vaf", "event"]
    )
    report_filename = "{}_report.tsv".format(sample_name)
    report_path = os.path.join(datahub.args.outdir, report_filename)
    result_df.to_csv(report_path, sep="\t", index=False)
    datahub.cleanup()


# for sv in VariantFile("/home/varelad/subset_709.vcf").fetch():
#     datahub = get_datahub(sv)
#     run(datahub)
#     print("DONE")

datahub = get_datahub()
run_option2(datahub)
