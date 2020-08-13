#!/usr/bin/env python

import argparse
import glob
import sys
import os
from os.path import abspath
# from my_functions import ImageSorter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", type=str, help=("List of files to scrape. Enclose any glob "
                                                 "patterns in quotes so they are passed unexpanded"))
    parser.add_argument("-o", metavar='outputPath', dest='outputPath', nargs='+', type=str,
                        help="Path to write the good/bad file lists too")
    # parser.add_argument("-j", metavar='joinKeys', dest='joinKeys', nargs='+', type=str,
    #                     help="Keys to return joined together.")
    # parser.add_argument("-l", metavar='libraryLocation', dest='libraryLocation', type=str,
    #                     help="Location of library for precomputed results.")
    # parser.add_argument("-p", action='store_true', default=False, dest='printPerFile',
    #                     help="Print all keys for each file?")
    # parser.add_argument("--noWarn", action='store_true', help="Suppress warnings for keys not in header?",
    #                     default=False, dest='noWarn')
    # parser.add_argument("--walk", action='store_true', help="Ignore path glob and walk whole tree for"
    #                     " all fits and fits.gz files",
    #                     default=False, dest='walk')
    # parser.add_argument("--oneFilePerDir", action='store_true', help="If walking, only take one file from"
    #                     " each directory",
    #                     default=False, dest='oneFilePerDir')

    args = parser.parse_args()
    files = args.files
    outputPath = args.outputPath[0]

    # if walk:
    #     for dirpath, dirnames, filenames in os.walk(args.files):
    #         for filename in [f for f in filenames if f.endswith(".fits") or f.endswith(".fits.gz")]:
    #             files.append(abspath(os.path.join(dirpath, filename)))
    #             if oneFilePerDir:
    #                 break
    # else:
    # files = [abspath(f) for f in glob.glob(args.files)]

    if not files:
        print('Found no files matching: ' + args.files, file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    outputGood = os.path.join(outputPath, 'good.txt')
    outputBad = os.path.join(outputPath, 'bad.txt')

    print(f"Output good: {outputGood}")
    print(f"Output bad: {outputBad}")

    # imageSorter = ImageSorter()

    # imageSorter.sort(files, outputGood, outputBad)


if __name__ == '__main__':
    main()
