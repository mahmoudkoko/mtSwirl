#!/usr/bin/env python3

import argparse

import hail as hl


def parse_args():
    parser = argparse.ArgumentParser(
        description="Lift homoplasmic variants onto the personalized mitochondrial reference."
    )
    parser.add_argument("--self-fasta", required=True)
    parser.add_argument("--self-fai", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--reference-fai", required=True)
    parser.add_argument("--reference-to-self-chain", required=True)
    parser.add_argument("--shifted-self-fasta", required=True)
    parser.add_argument("--shifted-self-fai", required=True)
    parser.add_argument("--forward-shift-chain", required=True)
    parser.add_argument("--filtered-vcf", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--output-filter-vcf", required=True)
    parser.add_argument("--output-shifted-vcf", required=True)
    return parser.parse_args()


def fai_to_len(fai):
    with open(fai, encoding="utf-8") as source:
        line = source.readline()
    return int(line.split("\t")[1])


def check_vcf_integrity(mt):
    if sorted(list(mt.row_key)) != ["alleles", "locus"]:
        raise ValueError("VCFs must always be keyed by locus, alleles.")

    if mt.aggregate_rows(~hl.agg.all(hl.len(mt.alleles) == 2)):
        raise ValueError(
            "This function only supports biallelic sites (run SplitMultiAllelics!)."
        )

    if mt.aggregate_rows(~hl.agg.all(hl.is_defined(mt.locus))):
        raise ValueError(
            "ERROR: locus must always be defined, both before and after Liftover. "
            "This should be a reversible operation, thus finding missing loci after "
            "reverse liftover is very concerning."
        )

    if mt.aggregate_rows(~hl.agg.all(hl.all(hl.map(hl.is_defined, mt.alleles)))):
        raise ValueError("ERROR: alleles should always be defined.")


def apply_conversion(mt, liftover_target, skip_flip=False):
    mt_lifted = mt.annotate_rows(
        new_locus=hl.liftover(mt.locus, liftover_target)
    )
    if skip_flip:
        mt_lifted = mt_lifted.annotate_rows(allele_flip=mt_lifted.alleles)
    else:
        mt_lifted = mt_lifted.annotate_rows(
            allele_flip=hl.reversed(mt_lifted.alleles)
        )
    mt_lifted = (
        mt_lifted.key_rows_by()
        .rename({"locus": "locus_orig", "alleles": "alleles_orig"})
        .rename({"new_locus": "locus", "allele_flip": "alleles"})
        .key_rows_by("locus", "alleles")
    )
    return mt_lifted.drop("locus_orig", "alleles_orig")


def main():
    args = parse_args()

    target = hl.ReferenceGenome(
        "target_self",
        ["chrM"],
        {"chrM": fai_to_len(args.self_fai)},
        mt_contigs=["chrM"],
    )
    source = hl.ReferenceGenome(
        "mtGRCh38",
        ["chrM"],
        {"chrM": fai_to_len(args.reference_fai)},
        mt_contigs=["chrM"],
    )
    target.add_sequence(args.self_fasta, args.self_fai)
    source.add_sequence(args.reference_fasta, args.reference_fai)
    source.add_liftover(args.reference_to_self_chain, "target_self")

    shifted_target = hl.ReferenceGenome(
        "target_self_shifted",
        ["chrM"],
        {"chrM": fai_to_len(args.shifted_self_fai)},
        mt_contigs=["chrM"],
    )
    shifted_target.add_sequence(args.shifted_self_fasta, args.shifted_self_fai)
    target.add_liftover(args.forward_shift_chain, shifted_target)

    mt_new = hl.import_vcf(args.filtered_vcf, reference_genome="mtGRCh38").select_entries()
    check_vcf_integrity(mt_new)

    mt_new_1 = mt_new.select_rows()
    mt_lifted_target = apply_conversion(mt_new_1, "target_self")
    check_vcf_integrity(mt_lifted_target)
    hl.export_vcf(mt_lifted_target, args.output_vcf, tabix=True)

    mt_new_2 = mt_new.select_rows("filters")
    this_metadata = hl.get_vcf_metadata(args.filtered_vcf)
    this_metadata = {"filter": this_metadata["filter"]}
    mt_lifted_target_withfilter = apply_conversion(mt_new_2, "target_self")
    check_vcf_integrity(mt_lifted_target_withfilter)
    hl.export_vcf(
        mt_lifted_target_withfilter,
        args.output_filter_vcf,
        tabix=True,
        metadata=this_metadata,
    )

    mt_lifted_shifted_target = apply_conversion(
        mt_lifted_target, "target_self_shifted", skip_flip=True
    )
    check_vcf_integrity(mt_lifted_shifted_target)
    hl.export_vcf(mt_lifted_shifted_target, args.output_shifted_vcf, tabix=True)


if __name__ == "__main__":
    main()
