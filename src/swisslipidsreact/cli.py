import logging
import logging.config
import argparse

from .utils import DEBUG
from .main import run_pipeline
from .ttl_export import export_ttl

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
    parser_run = subparsers.add_parser("run", help="Run the data processing pipeline")
    parser_run.add_argument(
        "--output-dir",
        type = str,
        default = None,
        help = "Output directory (default: current working directory)"
    )
    parser_run.add_argument(
        "--filter-fa",
        type = str,
        choices = ["curated", "c16", "none"],
        default = "curated",
        help = "Filter the fatty acids: curated (default), c16, none"
    )
    parser_run.add_argument(
        "--filter-rhea",
        type = str,
        choices = ["two-sides", "one-side"],
        default = "two-sides",
        help = "Filter Rhea by having an SLM class of an isomeric subspecies on at least one or both sides of the reaction: two-sides (default), one-side"
    )
    parser_run.add_argument(
        "--rhea-id",
        type = int,
        default = None,
        help = "Run pipeline only for the given Rhea ID"
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
        help = "Test run with palmitic acid only (default: False)"
    )
    
    # RDF export command
    parser_export = subparsers.add_parser("export-ttl", help = "Export RDF file from results")
    parser_export.add_argument(
        "--curated-fa",
        action = "store_true",
        help = "Use curated fatty acid list for TTL export (default: False for c16)"
    )
    parser_export.add_argument(
        "--input",
        type = str,
        default = None,
        help = "Input TSV file (default: inferred from mode)"
    )
    parser_export.add_argument(
        "--output-dir",
        type = str,
        default = None,
        help = "Output directory (default: current working directory)"
    )
    parser_export.add_argument(
        "--output-format",
        type = str,
        default = "ttl",
        help = "RDF output format for exporting (nt, ttl etc.)"
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
    elif args.command == "export-ttl":
        export_ttl(
            full_scope = args.curated_fa,
            input_path = args.input,
            output_dir = args.output_dir,
            output_format = args.output_format
        )
