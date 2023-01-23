#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reference kinematics associated functions and variables
@Author: Mads Jeppesen
@Date: 9/20/22
"""
import math
from scipy.spatial.transform import Rotation as R
from pyrosetta.rosetta.core.kinematics import Stub, Jump
from pyrosetta.rosetta.numeric import xyzMatrix_double_t, xyzVector_double_t
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num, jump_num_sym_dof
from pyrosetta import Pose
import random

dof_str_to_int = {
    "x" : 1,
    "y" : 2,
    "z" : 3,
    "angle_x" : 4,
    "angle_y" : 5,
    "angle_z" : 6,
    # FIXME: DELETE THIS IN FUTURE WHEN FIXED IN THE RIGIDBODYDOFADAPTIVEMOVER
    "x_angle": 4,
    "y_angle": 5,
    "z_angle": 6,
}

dof_int_to_str = {
    1: "x",
    2: "y",
    3: "z",
    4: "angle_x",
    5: "angle_y",
    6: "angle_z"
}

def get_x_rotation(jump: Jump):
    """Gets the x rotation associated with a Jump that only rotates around the x axis"""
    return R.from_matrix(jump.get_rotation()).as_euler("xyz", degrees=True)[0]

def get_y_rotation(jump: Jump):
    """Gets the y rotation associated with a Jump that only rotates around the y axis"""
    return R.from_matrix(jump.get_rotation()).as_euler("yxz", degrees=True)[0]

def get_z_rotation(jump: Jump):
    """Gets the z rotation associated with a Jump that only rotates around the z axis"""
    return R.from_matrix(jump.get_rotation()).as_euler("zyx", degrees=True)[0]

def get_x_translation(jump: Jump):
    """Gets the x translation associated with a Jump"""
    return jump.get_translation()[0]

def get_y_translation(jump: Jump):
    """Gets the y translation associated with a Jump"""
    return jump.get_translation()[1]

def get_z_translation(jump: Jump):
    """Gets the z translation associated with a Jump"""
    return jump.get_translation()[2]

def convert_into_rosetta(m):
    r = xyzMatrix_double_t(0)
    r.xx = m[0][0]
    r.xy = m[0][1]
    r.xz = m[0][2]
    r.yx = m[1][0]
    r.yy = m[1][1]
    r.yz = m[1][2]
    r.zx = m[2][0]
    r.zy = m[2][1]
    r.zz = m[2][2]
    return r

def x_rotation_matrix(angle, degrees=False):
    angle = angle * math.pi / 180 if degrees else angle
    return convert_into_rosetta([
        [1, 0, 0],
        [0, math.cos(angle), -math.sin(angle)],
        [0, math.sin(angle), math.cos(angle)]])

def y_rotation_matrix(angle, degrees=False):
    angle = angle * math.pi / 180 if degrees else angle
    return convert_into_rosetta([
        [math.cos(angle), 0, math.sin(angle)],
        [0, 1, 0],
        [-math.sin(angle), 0, math.cos(angle)]])

def z_rotation_matrix(angle, degrees=False):
    angle = angle * math.pi / 180 if degrees else angle
    return convert_into_rosetta([
        [math.cos(angle), -math.sin(angle), 0],
        [math.sin(angle), math.cos(angle), 0],
        [0, 0, 1]])

def set_all_dofs_to_zero(pose):
    """Sets all the movable dofs in the pose to zero"""
    for jump, dofs in get_dofs(pose).items():
        for dof, old_val in dofs.items():
            set_jumpdof_str_str(pose, jump, dof, 0)

# --- perturb from current position a particular dof value from jump.          --- #
# --- Below we have different function getters for that mainly for convenience --- #

def perturb_jumpdof_int_int(pose, jump: int, dof: int, value):
    assert is_symmetric(pose), "pose must be symmetrical"
    assert dof >= 1 and dof <= 6
    flexible_jump = pose.jump(jump)
    rot = flexible_jump.get_rotation()
    trans = flexible_jump.get_translation()
    if dof < 4:
        trans[dof - 1] += value
        flexible_jump.set_translation(trans)
    else:
        if dof == 4:
            rotmul = x_rotation_matrix(value, degrees=True)
        elif dof == 5:
            rotmul = y_rotation_matrix(value, degrees=True)
        else:
            rotmul = z_rotation_matrix(value, degrees=True)
        rot = rot * rotmul
        flexible_jump.set_rotation(rot)
    pose.set_jump(jump, flexible_jump)
    return pose

def perturb_jumpdof_str_int(pose, jump:str, dof : int, value: float):
    """Perturb the jump at dof of the symmetrical pose with a value"""
    jumpid = sym_dof_jump_num(pose, jump)
    return perturb_jumpdof_int_int(pose, jumpid, dof, value)

def perturb_jumpdof_str_str(pose, jump:str, dof : str, value: float):
    """Perturb the jump at dof of the symmetrical pose with a value"""
    jumpid = sym_dof_jump_num(pose, jump)
    return perturb_jumpdof_int_int(pose, jumpid, dof_str_to_int[dof], value)

# ------------------------------------------------ #

# --- Get a particular dof value from jump. Below we have different function getters for that mainly for convenience --- #

def get_jumpdof_from_jump(jump: Jump, dof: int):
    """Gets the value of a dof for a jump."""
    assert dof >= 1 and dof <= 6
    if dof == 1:
        return get_x_translation(jump)
    elif dof == 2:
        return get_y_translation(jump)
    elif dof == 3:
        return get_z_translation(jump)
    elif dof == 4:
        return get_x_rotation(jump)
    elif dof == 5:
        return get_y_rotation(jump)
    elif dof == 6:
        return get_z_rotation(jump)

def get_jumpdof_int_int(pose, jump: int, dof: int):
    """Gets the value of a dof using the jumpid directly (int) and dof (int)."""
    assert dof >= 1 and dof <= 6
    flexible_jump = pose.jump(jump)
    return get_jumpdof_from_jump(flexible_jump, dof)

def get_jumpdof_str_int(pose, jump: str, dof: int):
    """Gets the value of a dof using the str name of a jump (str) and the dof (int)."""
    assert dof >= 1 and dof <= 6
    jumpid = sym_dof_jump_num(pose, jump)
    return get_jumpdof_int_int(pose, jumpid, dof)

def get_jumpdof_str_str(pose, jump: str, dof: str):
    """Gets the value of a dof using the str name of a jump (str) and the dof (int)."""
    dof = dof_str_to_int[dof]
    assert dof >= 1 and dof <= 6
    jumpid = sym_dof_jump_num(pose, jump)
    return get_jumpdof_int_int(pose, jumpid, dof)

# ------------------------------------------------ #

# --- Set a particular dof value from jump. Below we have different function getters for that mainly for convenience --- #

def set_jumpdof_in_jump(jump: Jump, dof: int, value: float):
    """Set the dof in the jump with the exact value"""
    assert dof >= 1 and dof <= 6
    if dof < 4:
        trans = jump.get_translation()
        trans[dof - 1] = value
        jump.set_translation(trans)
    else:
        if dof == 4:
            rot = x_rotation_matrix(value, degrees=True)
        elif dof == 5:
            rot = y_rotation_matrix(value, degrees=True)
        else:
            rot = z_rotation_matrix(value, degrees=True)
        jump.set_rotation(rot)
    return jump

def set_jumpdof_int_int(pose, jump: int, dof : int, value: float):
    """Set the jump (with jumpid => int) at dof (int) of the symmetrical pose with the exact value"""
    assert is_symmetric(pose), "pose must be symmetrical"
    flexible_jump = pose.jump(jump)
    pose.set_jump(jump, set_jumpdof_in_jump(flexible_jump, dof, value))
    return pose

def set_jumpdof_str_int(pose, jump: str, dof : int, value: float):
    """Set the jump (with jumpname => str) at dof (int) of the symmetrical pose with the exact value"""
    assert is_symmetric(pose), "pose must be symmetrical"
    jump_id = sym_dof_jump_num(pose, jump)
    set_jumpdof_int_int(pose, jump_id, dof, value)
    return pose

def set_jumpdof_int_str(pose, jump: str, dof : int, value: float):
    """Set the jump (with jumpid => int) at dof (str) of the symmetrical pose with the exact value"""
    assert is_symmetric(pose), "pose must be symmetrical"
    set_jumpdof_int_int(pose, jump, dof_str_to_int[dof], value)
    return pose

def set_jumpdof_str_str(pose, jump: str, dof : str, value: float):
    """Set the jump (with jumpname => str) at dof (str) of the symmetrical pose with the exact value"""
    assert is_symmetric(pose), "pose must be symmetrical"
    jump_id = sym_dof_jump_num(pose, jump)
    set_jumpdof_int_int(pose, jump_id, dof_str_to_int[dof], value)
    return pose

def set_all_translations_to_0(pose):
    """Sets all translational dofs to 0"""
    for jump, params in get_dofs(pose).items():
        for dof, val in params.items():
            if not "angle" in dof:
                set_jumpdof_str_str(pose, jump, dof, 0)

def get_dofs(pose) -> dict:
    """Get the names of the movable jumps and their respective movable dof values as a dictionary.
    It look like this: { "jump name 1" : {"z" : 58.8}, "jump name 2" : {"y" : 42, "angle_x" : 32, .. } ..  }."""
    assert is_symmetric(pose), "pose must be symmetrical"
    position_info = {}
    dofs = pose.conformation().Symmetry_Info().get_dofs()
    for jump_id, symdof in dofs.items():
        jump_name = jump_num_sym_dof(pose, jump_id)
        flexible_jump = pose.jump(jump_id)
        for dof in range(1, 7):
            if symdof.allow_dof(dof):
                if position_info.get(jump_name) == None:
                    position_info[jump_name] = {}
                if dof == 1:
                    position_info[jump_name]["x"] = get_x_translation(flexible_jump)
                elif dof == 2:
                    position_info[jump_name]["y"] = get_y_translation(flexible_jump)
                elif dof == 3:
                    position_info[jump_name]["z"] = get_z_translation(flexible_jump)
                elif dof == 4:
                    position_info[jump_name]["angle_x"] = get_x_rotation(flexible_jump)
                elif dof == 5:
                    position_info[jump_name]["angle_y"] = get_y_rotation(flexible_jump)
                elif dof == 6:
                    position_info[jump_name]["angle_z"] = get_z_rotation(flexible_jump)
    return position_info

def get_dofs_as_values(pose) -> list:
    """Returns a list of the dofs in the pose as in the order they appear. See 'get_dofs' to retrieve it in a dictionary format."""
    assert is_symmetric(pose), "pose must be symmetrical"
    position_info = []
    dofs = pose.conformation().Symmetry_Info().get_dofs()
    for jump_id, symdof in dofs.items():
        jump_name = jump_num_sym_dof(pose, jump_id)
        flexible_jump = pose.jump(jump_id)
        for dof in range(1, 7):
            if symdof.allow_dof(dof):
                if dof == 1:
                    position_info.append(get_x_translation(flexible_jump))
                elif dof == 2:
                    position_info.append(get_y_translation(flexible_jump))
                elif dof == 3:
                    position_info.append(get_z_translation(flexible_jump))
                elif dof == 4:
                    position_info.append(get_x_rotation(flexible_jump))
                elif dof == 5:
                    position_info.append(get_y_rotation(flexible_jump))
                elif dof == 6:
                    position_info.append(get_z_rotation(flexible_jump))
    return position_info

def set_dofs(pose, jumpinfo):
    """Set the dofs from jumpdof info map such as the one that is outputted from 'get_dofs' above"""
    assert is_symmetric(pose), "pose must be symmetrical"
    for jump, params in jumpinfo.items():
        for doftype, value in params.items():
            set_jumpdof_str_str(pose, jump, doftype, value)


