#!/usr/bin/env python3

import sys, os, re


def main():

    usage = "\n\tusage: {} features_want_col1 feature.matrix\n\n".format(sys.argv[0])
    if len(sys.argv) < 3:
        exit(usage)

    features_want_file = sys.argv[1]
    matrix_file = sys.argv[2]

    features_want = set()
    with open(features_want_file) as fh:
        for line in fh:
            line = line.rstrip()
            feature_want = line.split("\t")[0]
            features_want.add(feature_want)

    with open(matrix_file) as fh:
        header = next(fh)
        print(header, end="")
        for line in fh:
            feature_id = line.split("\t")[0]
            if feature_id in features_want:
                print(line, end="")

    sys.exit(0)


if __name__ == "__main__":
    main()
