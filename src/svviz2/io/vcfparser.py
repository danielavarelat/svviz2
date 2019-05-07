import logging
import pysam
import re

from svviz2.utility.intervals import Locus
from svviz2.app import variants

logger = logging.getLogger(__name__)


class VCFParserError(Exception):
    pass


def only_nucs(seq):
    seq = seq.upper()
    return set(list(seq)) <= set(list("ACGT"))


def fix_vcf_header(vcf):
    if not "END" in vcf.header.info:
        # this is probably a bug in pysam, where it doesn't parse the END coordinate into variant.stop
        # if it's not defined in the header but doesn't let you read it through variant.info["END"]

        vcf.header.add_line(
            """##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate (exclusive)">"""
        )


class VCFParser(object):
    def __init__(self, datahub):
        self.datahub = datahub
        self.vcf = pysam.VariantFile(datahub.args.variants, drop_samples=True)
        fix_vcf_header(self.vcf)

    def get_variants(self):
        breakends = {}
        i = 0
        for variant in self.vcf:
            i = i + 1
            if not variant.id:
                variant.id = str(i)
            #                 raise VCFParserError("Variant ID must be specified in the VCF")

            if not "SVTYPE" in variant.info:
                if only_nucs(variant.ref) and only_nucs(
                    variant.alts[0]
                ):  # and sv_type == "INS":
                    yield get_sequence_defined(variant, self.datahub)
                    continue
                else:
                    print(
                        "variant does not appear to be a structural variant, skipping:{}".format(
                            variant
                        )
                    )
                    continue
            ######

            #             yield get_breakend(variant, variant, self.datahub)

            ####
            sv_type = variant.info["SVTYPE"].upper()
            if sv_type == "BND":
                mateid = variant.info["MATEID"]
                if not isinstance(mateid, str):
                    if len(mateid) > 1:
                        logger.error(
                            "ERROR: not sure what to do about this mateid: '{}'".format(
                                mateid
                            )
                        )
                    else:
                        mateid = mateid[0]

                if "," in mateid:
                    NotImplementedError(
                        "we currently don't support ambiguous breakends"
                    )

                if mateid in breakends:
                    breakend = get_breakend(
                        variant, breakends.pop(mateid), self.datahub
                    )
                    if breakend is not None:
                        yield breakend
                else:
                    assert not variant.id in breakends
                    breakends[variant.id] = variant
            elif only_nucs(variant.ref) and only_nucs(
                variant.alts[0]
            ):  # and sv_type == "INS":
                yield get_sequence_defined(variant, self.datahub)
            elif sv_type == "DEL":
                yield get_deletion(variant, self.datahub)
            elif sv_type == "INS" and (
                "INS:ME" in variant.alts[0] or "MEINFO" in variant.info
            ):
                raise NotImplementedError(
                    "not yet implemented: mobile element insertions"
                )
            elif sv_type == "TRA":
                yield get_translocation(variant, self.datahub)
            elif sv_type == "INV":
                yield get_inversion(variant, self.datahub)
            elif sv_type == "DUP":
                variant.alts == "'<DUP:TANDEM>': {}".format(variant)
                # if len(variant.alts) == 1 and variant.alts[0] == "<DUP:TANDEM>":
                yield get_tandem_duplication(variant, self.datahub)
                # else:
                #     logger.warn(
                #         "Only tandem duplications are supported; if this duplication is in fact "
                #         "a tandem duplication, make sure that the alt field of the vcf record "
                #         "is '<DUP:TANDEM>': {}".format(variant)
                #     )
            else:
                logger.warn("SKIPPING VARIANT: {}".format(variant))

        if len(breakends) > 0:
            logger.warn("found {} unpaired breakends".format(len(breakends)))


def get_sequence_defined(variant, datahub):
    # print("::", variant.id, variant.start, variant.stop, len(variant.ref))
    # variant.start: 0-based, inclusive; variant.stop: 0-based, exclusive
    if variant.stop - variant.start != len(variant.ref):
        error = (
            "VCF format error: coordinates ({}:{}-{}) do not match the variant length ({}). Please check the VCF variant"
            "spec; in particular, END coordinates are inclusive. Full variant: {}"
        )
        raise VCFParserError(
            error.format(
                variant.chrom, variant.start, variant.stop, variant.rlen, variant
            )
        )

    if len(variant.alts[0]) == 1 and variant.ref[0] == variant.alts[0][0]:
        # we need to add 1 to the start position to take into account the fact that the
        # alt is defined as the first nucleotide of the ref (which we've verified is actually
        # true here); we'll remove it from the ref so that the alt is zero-length
        deletion = variants.Deletion.from_breakpoints(
            variant.chrom, variant.start + 1, variant.stop - 1, datahub, variant.id
        )
        return deletion

    sdv = variants.SequenceDefinedVariant(
        variant.chrom,
        variant.start,
        variant.stop - 1,
        variant.alts[0],
        datahub,
        variant.id,
    )

    return sdv


def get_breakend(first, second, datahub):
    # return "{}\n{}".format(_parse_breakend(first), _parse_breakend(second))
    return parse_breakend(first, second, datahub)


def get_deletion(variant, datahub):
    # if variant.stop-1 > variant.start:
    #     stop = variant.stop
    # # elif "END" in variant.info:
    #     # stop = int(variant.info["END"])
    # else:
    #     errstr = "Error parsing event: '{}' -- missing 'END' coordinate; is END defined in the VCF header?"
    #     raise IOError(errstr.format("{}:{}-{} ({})".format(variant.chrom, variant.start, variant.stop, variant)))

    deletion = variants.Deletion.from_breakpoints(
        variant.chrom, variant.start, variant.stop - 1, datahub, variant.id
    )
    print("))))DEL:", deletion)
    return deletion


def get_tandem_duplication(variant, datahub):
    chrom, start, end = variant.chrom, variant.start - 1, variant.stop - 1
    strand = "+"
    duplicated_sequence = datahub.genome.get_seq(chrom, start, end, strand)

    sdv = variants.SequenceDefinedVariant(
        chrom, end, end, duplicated_sequence, datahub, variant.id
    )

    return sdv


def _parse_breakend(record):
    ref = record.ref
    alt = record.alts[0]
    if not ("[" in alt or "]" in alt):
        return None

    orientation = None
    altre1 = "(\[|\])(\w*):(\w*)(\[|\])(.*)"
    match = re.match(altre1, alt)
    if match:
        dir1, other_chrom, other_pos, dir2, alt_seq = match.groups()
        assert dir1 == dir2, alt
        if alt_seq != ref:
            raise Exception("not yet implemented: complex event")
        chrom, pos = record.chrom, record.pos

        if dir1 == "]":  # eg ]13:123456]T bnd_V
            orientation = "--"
        else:  # eg [17:198983[A bnd_X
            orientation = "-+"
    else:
        altre2 = "(.*)(\[|\])(\w*):(\w*)(\[|\])"
        match = re.match(altre2, alt)
        if match:
            alt_seq, dir1, other_chrom, other_pos, dir2 = match.groups()
            assert dir1 == dir2, (dir1, dir2)
            if alt_seq != ref:
                raise Exception("not yet implemented: complex event")
            chrom, pos = record.chrom, record.pos

            if dir1 == "]":  # eg G]17:198982 bnd_W
                orientation = "+-"
            else:  # eg C[2:321682[ bnd_U
                orientation = "++"

    if orientation is None:
        return None
    else:
        id_ = record.id
        if "EVENT" in record.info:
            # TODO: see if we care that a complex event can have multiple breakends
            # with the same "EVENT"
            id_ = record.info["EVENT"]

        result = {
            "chrom": chrom,
            "pos": pos,
            "other_chrom": other_chrom,
            "other_pos": int(other_pos),
            "orientation": orientation,
            "alt": alt,
            "id": id_,
        }
        return result


def parse_breakend(record1, record2, datahub):
    result1 = _parse_breakend(record1)

    result2 = _parse_breakend(record2)
    if not (
        result1["chrom"] == result2["other_chrom"]
        and result1["pos"] == result2["other_pos"]
    ):
        print(result1)
        print(result2)
        logger.error(
            "Malformed VCF: breakends do not appear to match:\n{}\n{}".format(
                record1, record2
            )
        )
        return None

    if (
        result1["chrom"] == result1["other_chrom"]
        and abs(result1["pos"] - result1["other_pos"]) < datahub.align_distance * 5
    ):
        logger.error("Can't yet handle nearby breakends; skipping")
        return None
    breakpoint1 = Locus(
        result1["chrom"],
        result1["pos"] - 1,
        result1["pos"] - 1,
        result1["orientation"][0],
    )
    breakpoint2 = Locus(
        result1["other_chrom"],
        result1["other_pos"] - 1,
        result1["other_pos"] - 1,
        result1["orientation"][1],
    )

    return variants.Breakend(breakpoint1, breakpoint2, datahub, result1["id"])


def get_inversion(record, datahub):
    return variants.Inversion(
        record.chrom, record.start, record.stop - 1, datahub, record.id
    )


def get_translocation(record, datahub):
    # breakpoint1 = Locus(record.chrom, record.start, record.start, "+")
    # breakpoint2 = Locus(
    #     record.info["CHR2"],
    #     record.stop,
    #     record.stop,
    #     "+" if record.info["STRANDS"][0] == "+" else "-",
    # )
    result1 = _parse_breakend(record)
    breakpoint1 = Locus(
        result1["chrom"],
        result1["pos"] - 1,
        result1["pos"] - 1,
        result1["orientation"][0],
    )
    breakpoint2 = Locus(
        result1["other_chrom"],
        result1["other_pos"] - 1,
        result1["other_pos"] - 1,
        result1["orientation"][1],
    )
    return variants.Breakend(breakpoint1, breakpoint2, datahub, record.id)


#     raise NotImplementedError()
