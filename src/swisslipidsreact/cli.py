import logging
import logging.config
import argparse

from .utils import DEBUG
from .main import run_pipeline
from .build_rdf import build_rdf

logger = logging.getLogger(__name__)

def main():

    # Configure logging:
    # * fixed width left-aligned, e.g. %(levelname)-8s
    # * truncate if too long, e.g. %(module)-20.20s
    #   (width = 20, precision = 20 -> truncate to 20 characters)
    logging.basicConfig(
        level = logging.DEBUG if DEBUG > 0 else logging.INFO,
        format = "%(asctime)s %(levelname)-8s - %(module)-20.20s %(lineno)4d: %(message)s",
        datefmt = "%Y-%m-%d %H:%M"
    )
    
    parser = argparse.ArgumentParser(description="SwissLipids Reaction Pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # default pipeline command
    parser_run = subparsers.add_parser("run", help = "Enumerate reactions")
    parser_run.add_argument(
        "--output-dir",
        type = str,
        default = None,
        help = "Output directory (default: current working directory)"
    )
    parser_run.add_argument(
        "--filter-fa",
        type = str,
        choices = ["c16", "curated", "none"],
        default = "c16",
        help = "Filter the fatty acids: c16 (default), curated, none (use only in combination with --rhea-id option)"
    )
    parser_run.add_argument(
        "--filter-rhea",
        type = str,
        choices = ["two-sides", "one-side"],
        default = "two-sides",
        help = "Filter Rhea by having a direct SLM parent class of an isomeric subspecies on at least one or both sides of the reaction: two-sides (default), one-side"
    )
    parser_run.add_argument(
        "--rhea-id",
        type = int,
        default = None,
        help = "Enumerate reactions only for the given Rhea ID"
    )
    parser_run.add_argument(
        "--rhea-version",
        type = int,
        default = None,
        help = "Use the given Rhea release version (default: latest release)"
    )
    parser_run.add_argument(
        "--test",
        action = "store_true",
        help = "Use only SwissLipids compounds whose FAs are all palmitate (default: False)"
    )
    
    # RDF build command.
    parser_build_rdf = subparsers.add_parser("build-rdf", help="Build RDF from reaction enumeration result file")
    parser_build_rdf.add_argument(
        "--input",
        type = str,
        default = None,
        help = "Input TSV file (default: <output-dir>/enumerated_reactions.tsv)"
    )
    parser_build_rdf.add_argument(
        "--output-dir",
        type = str,
        default = None,
        help = "Output directory (default: current working directory)"
    )
    parser_build_rdf.add_argument(
        "--output-format",
        type = str,
        default = "nt",
        help = "RDF serialization format (default: nt)"
    )

    args = parser.parse_args()

    if args.command == "run":
        run_pipeline(
            output_dir = args.output_dir,
            filter_fa = args.filter_fa,
            filter_rhea = args.filter_rhea,
            rhea_id = args.rhea_id,
            rhea_version = args.rhea_version,
            test = args.test
        )
    elif args.command == "build-rdf":
        build_rdf(
            input_file = args.input,
            output_dir = args.output_dir,
            output_format = args.output_format
        )
