#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility functions for the symmetryhandler
@Author: Mads Jeppesen
@Date: 5/13/21
"""
from symmetryhandler.symmetryhandler import SymmetrySetup
from pyrosetta.rosetta.core.pose.datacache import CacheableDataType
from shapedesign.src.utilities.pose import get_position_info
from io import StringIO
from pyrosetta.rosetta.core.pose import add_comment, dump_comment_pdb
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit
from pyrosetta import Pose
from pyrosetta.rosetta.core.conformation.symmetry import SymmData
from pyrosetta.rosetta.std import istringstream
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover

def symmetrize_pose_from_comment(pose):
    symmetry_file = pose.data().get_ptr(CacheableDataType.STRING_MAP).map()["SYMMETRY"].replace("|", "\n")
    s = SymmData()
    s.read_symmetry_data_from_stream(istringstream(symmetry_file))
    setup = SetupForSymmetryMover(s)
    setup.apply(pose)
    return pose

def get_updated_symmetry_file_from_pose(pose, initial_symmetry_file):
    pose = pose.clone() # we need this in order to not clear the datacache!
    setup = SymmetrySetup()
    setup.read_from_file(initial_symmetry_file)
    setup.update_from_pose(pose)
    return setup.make_symmetry_definition()

def get_symmetry_comment_from_pose(pose):
    return pose.data().get_ptr(CacheableDataType.STRING_MAP).map()["SYMMETRY"]

def get_symmetry_file_from_pose(pose):
    return get_symmetry_comment_from_pose(pose).replace("| ", "\n")

def create_symmetric_pose(pose, return_symmetry_file=False):
    """Creates a symmetric pose from the SYMMETRY comment in pdb file. Remember to initialize Rosetta
    with '-pdb_comments' to read the comment"""
    symmetry_file = StringIO(get_symmetry_file_from_pose(pose))
    s = SymmetrySetup()
    s.read_from_file(symmetry_file)
    s.make_symmetric_pose(pose)
    _ = symmetry_file.seek(0) # reuse file again
    add_symmetry_as_comment(pose, symmetry_file) # Adds to pose again because it is lost in the symmetrization
    if return_symmetry_file:
        return s

def set_all_dofs_to_zero(pose):
    """Sets all the movable dofs in the pose to zero"""
    from shapedesign.src.utilities.kinematics import set_jumpdof
    from shapedesign.src.utilities.pose import  dof_map
    for jump, dofs in get_position_info(pose, dictionary=True).items():
        for dof, old_val in dofs.items():
            set_jumpdof(pose, jump, dof_map[dof], 0)

def add_symmetry_as_comment(pose, symmetry_file):
    """Adds the symmetryfile as a comment to the pose."""
    if isinstance(symmetry_file, StringIO):
        f = symmetry_file
    else:
        f = open(symmetry_file, "r")
    # First create 2 parts of the symmetry file and fill in the set_dofs below
    symmetry_file_info = []
    string = ""
    for l in f.read().splitlines():
        if len(l) > 7 and l[0:7] == "set_dof":
            if not symmetry_file_info:  # empty
                symmetry_file_info.append(string + " | ")
                string = ""
        else:
            if string:
                string += f" | {l}"
            else:
                string += f"{l}"
    symmetry_file_info.append(string)
    string = symmetry_file_info[0]  # first part of the symmetry definition
    # set set_dof from the position_info
    for jumpname, dofs in get_position_info(pose, dictionary=True).items():
        string += f"set_dof {jumpname}"
        for dofname, dofval in dofs.items():
            string += f" {dofname}({dofval})"
        string += " | "
    string += symmetry_file_info[1]  # second part of the symmetry definition
    add_comment(pose, "SYMMETRY", string)

def update_symmetry_comment(pose):
    symmetry_file = StringIO(get_symmetry_file_from_pose(pose))
    pose.data().clear() #FIXME: removes the entire cache - do we want that?
    add_symmetry_as_comment(pose, symmetry_file)

def output_pose_with_symmetry_comment(pose, outname):
    "Outputs a pose with the current symmetry of the pose as a comment."
    update_symmetry_comment(pose)
    symmetry_comment = get_symmetry_comment_from_pose(pose)
    set_all_dofs_to_zero(pose)
    new_pose = Pose()
    extract_asymmetric_unit(pose, new_pose, False)
    add_comment(new_pose, "SYMMETRY", symmetry_comment)
    dump_comment_pdb(outname, new_pose)
