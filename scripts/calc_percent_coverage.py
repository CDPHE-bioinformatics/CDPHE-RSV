#! /usr/bin/env python

"""Calculate percent coverage."""

import argparse
import logging
import numpy as np
import pandas as pd
import re
import sys

from Bio import SeqIO

__author__ = "CDPHE"
__copyright__ = "State of Colorado"
__license__ = "GPL-3.0.0-or-later"

log = logging.getLogger(__name__)


def parse_args(args) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate percent coverage.")

    parser.add_argument(
        "--log_level",
        help="the level to log at",
        choices=["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"],
        default="INFO",
    )
    parser.add_argument("--sample_name", help="sample name")
    parser.add_argument("--fasta_file", help="fasta file")
    parser.add_argument("--reference_file", help="reference file")

    return parser.parse_args(args)


def setup_logging(log_level) -> None:
    """Set up logging."""
    log_format = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=log_level,
        stream=sys.stdout,
        format=log_format,
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def calculate_percent_coverage(sample_name, fasta_file, reference_file):
    """Calculate percent coverage."""
    # first check that there is only one sequence in the fasta files
    num_records = 0
    with open(fasta_file, "r") as fasta_handle:
        inside_text = fasta_handle.read()
    for _ in re.finditer(">", inside_text):
        num_records = num_records + 1

    # calculate the percent coverage
    ref_fasta = SeqIO.read(reference_file, "fasta")
    ref_fasta_length = len(ref_fasta.seq)

    Ns = 0
    coverage = 0.0
    aligned_bases = 0
    num_non_ambiguous_bases = 0
    if num_records == 1:
        record = SeqIO.read(fasta_file, "fasta")
        aligned_bases = len(record.seq)
        if aligned_bases == 0:
            Ns = ref_fasta_length
            num_non_ambiguous_bases = 0
            coverage = 0
        else:
            uncorrected_Ns = record.seq.count("N")
            Ns = (ref_fasta_length - aligned_bases) + uncorrected_Ns
            num_non_ambiguous_bases = aligned_bases - uncorrected_Ns
            coverage = round((1 - (Ns / ref_fasta_length)) * 100, 2)
    else:
        aligned_bases = np.NaN  # type: ignore
        Ns = np.NaN  # type: ignore
        num_non_ambiguous_bases = np.NaN  # type: ignore
        coverage = np.NaN  # type: ignore

    # create pd df with calc_percent_cvg
    df = pd.DataFrame()
    df["sample_name"] = [sample_name]
    df["aligned_bases"] = [aligned_bases]
    df["N_bases"] = [Ns]
    df["non_ambiguous_bases"] = [num_non_ambiguous_bases]
    df["percent_coverage"] = [coverage]
    # df['number_seqs_in_fasta'] = [num_records]

    outfile = "%s_consensus_cvg_stats.csv" % sample_name
    df.to_csv(outfile, index=False)


def main(args):
    """Main function."""
    setup_logging(log_level=args.log_level)

    calculate_percent_coverage(
        sample_name=args.sample_name,
        fasta_file=args.fasta_file,
        reference_file=args.reference_file,
    )


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args=args)
