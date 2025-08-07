import argparse
from swisslipidsreact.main import run_pipeline
from swisslipidsreact.ttl_export import export_ttl
from swisslipidsreact.MasterIdAnalysis import MasterIdAnalysis

def main():
    parser = argparse.ArgumentParser(description="SwissLipids Reaction Pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # default pipeline command
    parser_run = subparsers.add_parser("run", help="Run the data processing pipeline")
    parser_run.add_argument(
        "--curated-fa",
        action="store_true",
        help="Use curated fatty acid list (default: False for C16)"
    )
    parser_run.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: current working directory)"
    )
    parser_run.add_argument(
        "--rheaid",
        type=int,
        default=None,
        help="run pipeline for only one rhea id"
    )
    parser_run.add_argument(
        "--all-fa",
        action="store_true",
        default=False,
        help="No restrictions of FA per position"
    )
    # ttl export command
    parser_export = subparsers.add_parser("export-ttl", help="Export RDF Turtle file from results")
    parser_export.add_argument(
        "--curated-fa",
        action="store_true",
        help="Use curated fatty acid list for TTL export (default: False for C16)"
    )

    parser_export.add_argument(
        "--input",
        type=str,
        default=None,
        help="Input TSV file (default: inferred from mode)"
    )
    parser_export.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: current working directory)"
    )

    parser_master_id_analysis = subparsers.add_parser("master-id-analysis", help="Export RDF Turtle file from results")
    parser_master_id_analysis.add_argument(
        "--input",
        type=str,
        default=None,
        help="Input TSV file (default: inferred from mode)"
    )
    parser_master_id_analysis.add_argument(
        "--curated-fa",
        action="store_true",
        help="Use curated fatty acid list (default: False for C16)"
    )
    parser_master_id_analysis.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: current working directory)"
    )
    parser_master_id_analysis.add_argument(
        "--all-fa",
        action="store_true",
        default=False,
        help="No restrictions of FA per position"
    )
    args = parser.parse_args()

    if args.command == "run":
        run_pipeline(
            curated_fa_list_run=args.curated_fa,
            output_dir=args.output_dir,
            no_curated_list_restrictions=args.all_fa,
            rheaid=args.rheaid
        )
    elif args.command == "export-ttl":
        export_ttl(
            full_scope=args.curated_fa,
            input_path=args.input,
            output_dir=args.output_dir
        )

    elif args.command == "master-id-analysis":
        analysis = MasterIdAnalysis(
            output_dir=args.output_dir
        )
        analysis.run_master_id_analysis(
            results_overview_path=args.input,
            curated_fa_list_run=args.curated_fa,
            output_dir=args.output_dir,
            no_curated_list_restrictions=args.all_fa
            )