#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_symmetry_file.py
@Author: Mads Jeppesen
@Date: 25/6/2024
"""
from pyrosetta import init
from symmetryhandler.symmetrymaker import SymmetryMaker

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Creates a Rosetta reference symmetry definition file for any point group symmetry")
    parser.add_argument(dest='struct', help="input file to symmetrize", type=str)
    parser.add_argument(dest='point_group', help="point group to make", type=str)
    parser.add_argument(dest='out', help="Directory path to output the symmetry files to", type=str)
    parser.add_argument('-chains', help="Chains to include in the symmetry file. Either a minimal "
                                        "representation (minimal) or a full representation (full).", type=str,
                                        choices=("minimal", "unique", "full"), default="full")
    parser.add_argument('-a', help="make_symmdef_file commands", nargs="+", type=str)
    parser.add_argument('-i', help="make_symmdef_file commands", nargs="+", type=str)
    parser.add_argument('-b', help="make_symmdef_file commands", nargs="+", type=str)
    args = parser.parse_args()

    init()

    sm = SymmetryMaker(args.struct, args.point_group, args.chains, args.out, args.a, args.i, args.b)
    if args.point_group == "cyclical":
        sm.make_cyclical_symmetry()

if __name__ == "__main__":
    main()