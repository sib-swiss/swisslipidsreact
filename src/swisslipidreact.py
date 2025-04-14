#!/usr/bin/env python3

import argparse
from swisslipidreact import RheaToSwisslipidsDf

def main():
    parser = argparse.ArgumentParser(description="CLI tool for RheaDB reaction utilities")
    parser.add_argument('--enumerate', action='store_true', help='')
    parser.add_argument('--SLMS', type=str, help='list of SLMs of the fatty acids, comma separated')

    args = parser.parse_args()

    r2sl_df = RheaToSwisslipidsDf()

    if args.enumerate:
        r2sl_df
    else:
        parser.print_help()

if __name__ == '__main__':
    main()