#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Converts a denovo symmetry script to a regular symmetry script
@Author: Mads Jeppesen
@Date: 6/17/22
"""

import argparse
from symmetryhandler.symmetryhandler import SymmetrySetup, CoordinateFrame
from symmetryhandler.mathfunctions import rotation_matrix
from copy import deepcopy
from pathlib import Path

def main(symmetry_file, outname=None):
    ss = SymmetrySetup()
    read_transformations = False
    with open(args.symmetry_file, "r") as f:
        for line in f:
            line = line.split()
            if read_transformations:
                if line[0] == "start":
                    x, y, o = [[int(i) for i in ii.split(',')] for ii in line[1:]]
                    ss.add_vrt(CoordinateFrame("VRT0001", x, y, orig=o, generate_z_from_xy=True))
                elif line[0] == "rot":
                    rot_type = line[1]
                    rot_times = int(line[2])
                    rot_degrees = 360 / rot_times
                    if rot_type == "Rz":
                        axis = ss.get_vrt_name("VRT0001").vrt_z
                    elif rot_type == "Ry":
                        axis = ss.get_vrt_name("VRT0001").vrt_y
                    elif rot_type == "Rx":
                        axis = ss.get_vrt_name("VRT0001").vrt_x
                    else:
                        raise NotImplementedError("Only Rz, Ry and Rx has code for it! Sorry")
                    for n, r in enumerate(range(1, rot_times), 2): # we have already done the first (0), and we are naming from 2!
                        new_vrt = deepcopy(ss.get_vrt_name("VRT0001"))
                        new_vrt.name = f"VRT{str(n).rjust(4, '0')}"
                        new_vrt.rotate(rotation_matrix(axis, rot_degrees * r))
                        ss.add_vrt(new_vrt)
            if line[0] == "symmetry_name":
                ss.symmetry_name = " ".join(line[1:])
            elif line[0] == "E":
                ss.energies = " ".join(line[2:])
            elif line[0] == "anchor_residue":
                ss.anchor = " ".join(line[1:])
            elif line[0] == "recenter":
                ss.recenter = True
            # specific for denovo
            elif line[0] == "virtual_transforms_start":
                read_transformations = True
                # make a global vrt
                ss.add_vrt(CoordinateFrame("VRTglobal", [1,0,0], [0,1,0], orig=[0, 0, 0]))
            elif line[0] == "virtual_transforms_stop":
                read_transformations = False
            elif line[0] == "set_dof":
                # we need to connect all vrts to global one
                for n, vrt in enumerate(ss._vrts[1:], 1): # skipt VRTglobal
                    ss.add_jump(f"JUMPG{n}", "VRTglobal", vrt.name)
                # then we connect the other vrts to the subunits
                for n, vrt in enumerate(ss._vrts[1:], 1): # skipt VRTglobal
                    ss.add_jump(f"JUMP{n}", vrt.name, "SUBUNIT")
                # now get the available dofs
                for dof in line[2:]:
                    if dof[:5] == "angle":
                        degree = "rotation"
                        axes = dof[6]
                    else:
                        degree = "translation"
                        axes = dof[0]
                    value = dof.split("(")[1].split(")")[0]
                    if ":" in value:
                        print("':' is not understood and will be left out")
                        value = 0
                    ss.add_dof("JUMP1", axes, degree, float(value))
                ss.add_jumpgroup("MODIFIED_BASEJUMP", *[j for j in ss._jumps if not "G" in j])
    if not outname:
        outname = Path(symmetry_file).stem + f"_converted.{Path(symmetry_file).suffix}"
    else:
        ss.output(outname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts a denovo symmetry script to a regular symmetry script")
    parser.add_argument('--symmetry_file', help="denovo symmetry file (contains the 'virtual_transforms_start' token)",
                        type=str, required=True)
    parser.add_argument('--outname', help="the name of the resulting outfile", type=str)
    args = parser.parse_args()
    main(args.symmetry_file, args.outname)

# symmetry_name c4
# subunits 4
# recenter
# number_of_interfaces  2
# E = 4*VRT0001 + 4*(VRT0001:VRT0002) + 2*(VRT0001:VRT0003)
# anchor_residue COM
# virtual_transforms_start
# start -1,0,0 0,1,0 0,0,0
# rot Rz 4
# virtual_transforms_stop
# connect_virtual JUMP1 VRT0001 VRT0002
# connect_virtual JUMP2 VRT0002 VRT0003
# connect_virtual JUMP3 VRT0003 VRT0004
# set_dof BASEJUMP x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)

# symmetry_name c4
# E = 4*VRT0001 + 4*(VRT0001:VRT0002) + 2*(VRT0001:VRT0003)
# anchor_residue COM
# recenter
# virtual_coordinates_start
# VRTglobal ...
# VRT0001
# VRT0002
# VRT0003
# VRT0004
# virtual_coordinates_stop
# connect_virtual JUMPglobal1 VRTglobal VRT0001
# connect_virtual JUMPglobal2 VRTglobal VRT0002
# connect_virtual JUMPglobal3 VRTglobal VRT0003
# connect_virtual JUMPglobal4 VRTglobal VRT0004
# set_dof JUMPglobal1 x(50) angle_x(0:360) angle_y(0:360) angle_z(0:360)
# set_jump_group JUMPglobal1 JUMPglobal2 JUMPglobal3 JUMPglobal4