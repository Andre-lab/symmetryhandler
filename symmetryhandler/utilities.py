#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility functions.
@Author: Mads Jeppesen
@Date: 5/13/21
"""
from symmetryhandler.symmetrysetup import SymmetrySetup
from pyrosetta.rosetta.core.pose.datacache import CacheableDataType
from io import StringIO
from pyrosetta.rosetta.core.pose import add_comment, dump_comment_pdb
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit
from pyrosetta import Pose
from pyrosetta.rosetta.core.conformation.symmetry import SymmData
from pyrosetta.rosetta.std import istringstream
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
# from symmetryhandler.kinematics import set_jumpdof, dof_map

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
    setup.update_dofs_from_pose(pose)
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


