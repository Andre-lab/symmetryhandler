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
    parser.add_argument(dest='point_group', help="point group to make", type=str, choices=("C", ))
    parser.add_argument(dest='out', help="Directory path to output the symmetry files to", type=str)
    parser.add_argument('-chains', help="Chains to include in the symmetry file. Either a minimal "
                                        "representation (minimal), unique (unique) or a full representation (full).", type=str,
                                        choices=("minimal", "unique", "full"), default="full")
    parser.add_argument('-align_struct', help="Align the structure to this pdb. Example 'file.pdb-A' will align "
                                              "the output symmetry to file.pdb and the main chain onto A.", type=str)
    parser.add_argument('-custom_energy', help="Custom energy line", type=str)
    parser.add_argument('-a', help="make_symmdef_file.pl commands", nargs="+", type=str)
    parser.add_argument('-i', help="make_symmdef_file.pl commands", nargs="+", type=str)
    parser.add_argument('-b', help="make_symmdef_file.pl commands", nargs="+", type=str)
    args = parser.parse_args()

    init()

    sm = SymmetryMaker(args.struct, args.point_group, args.chains, args.out, args.a, args.i,
                       args.b, args.align_struct, args.custom_energy)
    if args.point_group.upper() == "C":
        sm.make_cyclical_symmetry()

if __name__ == "__main__":
    main()