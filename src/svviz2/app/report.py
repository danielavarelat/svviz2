import collections
import json
import numpy
import os
import pandas

from svviz2.app.variants import non_negative
from svviz2.utility import statistics
from svviz2.remap import genotyping


def report(datahub):
    results = []

    results.extend(tally_support(datahub))
    # results.extend(tally_segments(datahub))
    # results.extend(tally_nearby_polymorphisms(datahub))

    result_df = pandas.DataFrame(results, columns=["sample", "allele", "key", "value"])

    result_df["event"] = datahub.variant.name
    report_filename = "{}.report.tsv".format(datahub.variant.short_name())
    report_path = os.path.join(datahub.args.outdir, report_filename)
    result_df.to_csv(report_path, sep="\t", index=False)


def tally_support(datahub):
    results = []
    for sample_name, sample in datahub.samples.items():
        allele_count = {}
        allele_mapq_sum = {}
        # allele_weighted_count = {}
        for allele in ["alt", "ref"]:
            bam = sample.outbam(allele, "r")
            cur_results, count, mapq_sum, _ = _tally_support(bam)
            allele_count[allele] = count
            allele_mapq_sum[allele] = mapq_sum
            # allele_weighted_count[allele] = weighted_count
            for cur_result in cur_results:
                results.append((sample_name, allele) + cur_result)

        # for __, value in [
        #     ("count", allele_count),
        # ("mapq", allele_mapq_sum),
        # ("weighted", allele_weighted_count),
        # ]:
        # print(name, value)
        # cur_genotype = genotyping.calculate_genotype_likelihoods(
        #     value["ref"], value["alt"]
        # )
        # results.append((sample_name, "", "GL_{}".format(name), cur_genotype[0]))
        # results.append((sample_name, "", "GQ_{}".format(name), cur_genotype[1]))

        # gt = max(zip(cur_genotype[0], ["0/0", "0/1", "1/1"]), key=lambda x: x[0])[1]
        # results.append((sample_name, "", "GT_{}".format(name), gt))
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
    # results = [("count", count), ("weighted_count", weighted_count)]

    # breakpoint_ids = {}
    # for i, breakpoint in enumerate(sorted(breakpoints)):
    #     breakpoint_ids[breakpoint] = i

    # for breakpoint, overlaps in breakpoint_overlaps.items():
    #     key = "overlap_{}.{}".format(breakpoint_ids[breakpoint], breakpoint)
    #     results.append((key, numpy.mean(overlaps)))

    # for breakpoint, cur_count in breakpoint_counts.items():
    #     breakpoint, overlaps_sequence = breakpoint
    #     key = "count_{}.{}_{}".format(
    #         breakpoint_ids[breakpoint],
    #         breakpoint,
    #         "seq" if overlaps_sequence else "pair",
    #     )
    #     results.append((key, cur_count))

    # for breakpoint, cur_extensions in extensions.items():
    #     key = "extension_{}.{}".format(breakpoint_ids[breakpoint], breakpoint)
    #     results.append((key, numpy.mean(cur_extensions)))

    return results, count, mapq_sum, weighted_count

