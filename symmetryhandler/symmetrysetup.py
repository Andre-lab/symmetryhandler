#!/usr/bin/env python3
# coding=utf-8
"""
SymmetrySetup class
@Author: Mads Jeppesen
@Date: 9/21/22
"""
import copy
import textwrap
import numpy as np
import xmlrpc.client as xmlrpclib
from pyrosetta.rosetta.std import istringstream
from pyrosetta.rosetta.core.conformation.symmetry import SymmData
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.pose.datacache import CacheableDataType
from pyrosetta import Pose
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit
from pyrosetta.rosetta.core.conformation.symmetry import residue_center_of_mass
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.core.kinematics import Stub, Jump
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num, jump_num_sym_dof
from symmetryhandler.coordinateframe import CoordinateFrame
from symmetryhandler.mathfunctions import rotation_matrix, rotate
from symmetryhandler.reference_kinematics import set_jumpdof_str_str, get_dofs
from scipy.spatial.transform import Rotation as R
from pathlib import Path
from io import StringIO
import yaml
from pyrosetta.rosetta.std import istringstream

class SymmetrySetup:
    """A symmetric setup that stores the symmetry information used internally in Rosetta.

    See this link to the Rosetta user page for a description of symmetry.:
    https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry#symmetry-definitions

    Has additional methods for what I call "reference symmetry".
    That is, a symmetry definition where the following conditions apply:
        1. Each degree of freedom has 1 dedicated and unique jump.
        2. Each jump have their 2 associated VRTs (stubs) in the same coordinate frame.
           That is, the x, y, z directions are the same but not nescesarily the origo.
    With such a definition the changes in the degrees of freedom can be tracked.
    """

    def __init__(self, file=None, pose=None, symmetry_name=None):
        """Initializes the Symmetric setup.

        :param str symmetry_name: name of the symmetric setup.
        :param str anchor: the residue to attach anchor the symmetry system to.
        :param str energies: the energies contained in the symmetric system.
        :param bool recenter: if to recenter the input protein before applying symmetry.
        :param list vrts: the list of CoordinateFrames in the symmetric system.
        :param dict jumps: the jumps in the in the symmetric system.
        :param dict jumpgroups: the jumpgroups in the symmetric system.
        :param dict dofs: the degrees of freedom in the symmetric system.
        """
        self.symmetry_name = symmetry_name
        self.anchor = None
        self.recenter = None
        self.energies = []
        self._jumps = {}
        self._jumpgroups = {}
        self._dofs = {}
        self._dofrefs = {}
        self._vrts = []
        self._init_vrts = None
        self.reference_symmetric = None
        self.actual_anchor_residue = None
        self.headers = {}
        if file:
            self.read_from_file(file)
        elif pose:
            self.read_from_pose_comment(pose)

    def get_map_pose_resi_to_vrt(self, pose):
        return {v: k for k, v in self.get_map_vrt_to_pose_resi(pose).items()}

    def get_map_vrt_to_pose_resi(self, pose):
        """Maps the virtual residue names to the pose residue numbers."""
        vrt_pose_resi = {}
        for jump, (vrt_up, vrt_down) in self._jumps.items():
            vrt_up_resi = pose.fold_tree().upstream_jump_residue(sym_dof_jump_num(pose, jump))
            vrt_down_resi = pose.fold_tree().downstream_jump_residue(sym_dof_jump_num(pose, jump))
            if not vrt_up in vrt_pose_resi:
                vrt_pose_resi[vrt_up] = vrt_up_resi
            if not vrt_down in vrt_pose_resi and not vrt_down == "SUBUNIT":
                vrt_pose_resi[vrt_down] = vrt_down_resi
        return vrt_pose_resi


    def vrts_overlap_with_pose(self, pose, update_and_apply_dofs=True, atol=1e-1):
        if update_and_apply_dofs:
            self.update_dofs_from_pose(pose, apply_dofs=True, assert_correct_overlap=False)
        for vrtname, resi in self.get_map_vrt_to_pose_resi(pose).items():
            residue = pose.residue(resi)
            assert residue.name() == "VRT"
            # get the residue coordinates
            r_o = np.array(residue.atom(1).xyz())
            r_x = np.array(residue.atom(2).xyz())
            r_y = np.array(residue.atom(3).xyz())
            # check they overlap with the one in the SymmetrySetup
            vrt = self.get_vrt(vrtname)
            v_o = vrt._vrt_orig
            v_x = vrt._vrt_x + vrt._vrt_orig
            v_y = vrt._vrt_y + vrt._vrt_orig
            try:
                assert np.isclose(r_o, v_o, atol=atol).all()
                assert np.isclose(r_x, v_x, atol=atol).all()
                assert np.isclose(r_y, v_y, atol=atol).all()
            except AssertionError:
                raise AssertionError
        return True


    @staticmethod
    def make_asymmetric_pose(pose_in, reset_dofs=True, dont_reset: list = None):
        """Extact the asymmetric pose with all reset dofs to 0."""
        assert is_symmetric(pose_in), "Pose is not symmetric!"
        pose = pose_in.clone()
        if reset_dofs:
            # assert self.reference_symmetric, ""
            # set all degrees of freedom to 0
            for jump, dofs in get_dofs(pose).items():
                if dont_reset and jump in dont_reset:
                    continue
                for dof, old_val in dofs.items():
                    set_jumpdof_str_str(pose, jump, dof, 0)
            # check that the anchor atom is actually at zero!
            com_ca = np.array(pose.residue(residue_center_of_mass(pose.conformation(), 1, pose.chain_end(1))).atom("CA").xyz())
            try:
                assert np.isclose(com_ca, [0, 0, 0], atol=1e-2).all()
            except:
                raise ValueError
        # create a new pose object and fill it with the asymmetric pose
        apose = Pose()
        extract_asymmetric_unit(pose, apose, False)
        return apose

    def get_trans_from_jump(self, jump, dofname):
        if dofname == "x":
            return jump.get_translation()[0]
        elif dofname == "y":
            return jump.get_translation()[1]
        elif dofname == "z":
            return jump.get_translation()[2]
        else:
            raise ValueError("Only dofname = 'x', 'y' or 'z' are understood.")

    def get_rot_from_jump(self, jump, dofname):
        if dofname == "angle_x":
            return R.from_matrix(jump.get_rotation()).as_euler("xyz", degrees=True)[0]
        elif dofname == "angle_y":
            return R.from_matrix(jump.get_rotation()).as_euler("yxz", degrees=True)[0]
        elif dofname == "angle_z":
            return R.from_matrix(jump.get_rotation()).as_euler("zyx", degrees=True)[0]
        else:
            raise ValueError("Only dofname == 'angle_x', 'angle_y' or 'angle_z' are understood.")

    def is_reference_symmetry(self):
        """Checks if the symmetry file is reference based."""
        is_reference_symmetry = True
        # 1. condition: check that all dofs have unique jumps
        jumps_used = []
        for jumpname in self._dofs.keys():
            if jumpname in jumps_used:
                is_reference_symmetry = False
            else:
                jumps_used.append(jumpname)
        # 2. condition: check that the vrts have the same coordinate frame
        for jumpname, params in self._dofs.items():
            down, up = self._jumps[jumpname]
            if up == "SUBUNIT":
                continue
            downstream_vrt = self.get_vrt(down)
            upstream_vrt = self.get_vrt(up)
            if params[0][1] == "rotation":
                if params[0][0] == "x":
                    is_reference_symmetry *= upstream_vrt.identical_x(downstream_vrt)
                elif params[0][0] == "y":
                    is_reference_symmetry *= upstream_vrt.identical_y(downstream_vrt)
                elif params[0][0] == "z":
                    is_reference_symmetry *= upstream_vrt.identical_z(downstream_vrt)
            else:  # translation
                if params[0][0] == "x":
                    is_reference_symmetry *= upstream_vrt.identical_x(downstream_vrt)
                elif params[0][0] == "z":
                    is_reference_symmetry *= upstream_vrt.identical_z(downstream_vrt)
        return is_reference_symmetry

    def reset_all_dofs(self):
        """Sets all the values of every dof to 0."""
        for jn, vl in self._dofs.items():
            for v in vl:
                v[2] = 0.0

    def reset_jumpdofs(self, jumpname):
        """Sets all the values of dofs belonging to a Jump to 0."""
        for dof in self._dofs[jumpname]:
            dof[2] = 0.0

    def set_dof(self, jumpname, dof, doftype, value):
        """Sets the value of the dof.

        :param jumpname: Jumpname
        :param dof: dof (x,y,z)
        :param doftype: (translation, rotation)
        :param value: value to set for the dof
        :return:
        """
        assert "angle" not in dof, "Only use 'x', 'y' or 'z'."
        for jn, vl in self._dofs.items():
            if jn == jumpname:
                for n, v in enumerate(vl):
                    t_dof, t_doftype, _ = v
                    if t_dof == dof and t_doftype == doftype:
                        self._dofs[jn][n][2] = value
                        return
        raise ValueError(f"Cannot set {jumpname}, {dof}, {doftype} because it does not exist!")

    def add_vrt(self, vrt):
        """Adds a CoordinateFrame to this instance.

        :param symmetryhandler.coordinateframe.CoordinateFrame vrt: name of the CoordinateFrame.
        """
        self._vrts.append(vrt)

    def remove_vrt(self, name):
        """Removes virtual residue coordinate systems (vrts) from the symmetry setup object.

        :param str name: name of vrt.
        """
        for vrt in self._vrts:
            if name == vrt:
                self._vrts.remove(name)

    def add_jump(self, name, src_name, dest_name):
        """Creates a jump between shapedock (source) and dest (destination).

        vrts must already be added with the add_vrt function, except for "SUBUNIT".

        :param str name: name of jump.
        :param str src_name: name of source vrt.
        :param str dest_name: name of destination vrt.
        """
        if name in self._jumps:
            raise ValueError(f"Name of jump {name} does already exist")
        vrt_names = [vrt.name for vrt in self._vrts] + ['SUBUNIT']
        if not src_name in vrt_names:
            raise ValueError(f"Name of source vrt: {src_name} does not exist")
        if not dest_name in vrt_names:
            raise ValueError(f"Name of destination vrt: {dest_name} does not exist")
        self._jumps[name] = [src_name, dest_name]

    def remove_jump(self, name):
        """Removes a jump in the tree with a given name.

        Note: This can create a disconnection in the tree.

        :param str name: name of the jump.
        """
        if not name in self._jumps:
            raise ValueError(f"Name of jump: {name} does not exist")
        self._jumps.pop(name)

    def add_jumpgroup(self, name, *jump_names):
        """Creates a jumpgroup.

        jumps must already have been added with the add_jumps function.

        :param str name: name of jumpgroup.
        :param str jump_names: names of jumps to include in the jumpgroup. The first name will be the master jump.
        """
        if name in self._jumpgroups:
            raise ValueError(f"Name of jumpgroup {name} already exist")
        jumps = self._jumps.keys()
        for jump_name in jump_names:
            if jump_name not in jump_names:
                raise ValueError("\"{}\" does not exist".format(jump_name))
        self._jumpgroups[name] = list(jump_names)

    def remove_jumpgroup(self, name):
        """Removes a jumpgroup with a given name.

        :param str name: name of the jumpgroup.
        """
        if not name in self._jumpgroups:
            raise ValueError(f"Name of jump: {name} does not exist")
        self._jumpgroups.pop(name)

    def add_dof(self, name, axes, degree, value=None):
        """Creates a degree of freedom for a jump.

        Jump must already be added with the add_vrt function.

        :param str name: name of the jump.
        :param str axes: the axes to apply the jump to. Must be either 'x', 'y' and 'z'.
        :param str degree: the degree of freedom in question. It's either 'translation' or 'rotation'.
        :param float value: the value of the degree of freedom.
        Ångström is assumed for translation and degrees for rotation.
        """
        if not name in self._jumps:
            raise ValueError(f"Name of jump: {name} does not exist")

        if not axes in 'xyz':
            raise ValueError("axes should either be 'x', 'y' or 'z' ")
        if not degree in ("translation", "rotation"):
            raise ValueError("axes should either be \"translation\", \"rotation\" ")
        if not name in self._dofs.keys():
            self._dofs[name] = []
        self._dofs[name].append([axes, degree, value])

    def remove_dofs(self, name):
        """Removes dofs for a given jump.

        :param str name: name of the jump
        """
        if not name in self._dofs:
            raise ValueError(f"Name of jump: {name} does not exist")
        self._dofs.pop(name)

    def get_vrt(self, name):
        """Returns the CoordinateFrame with the given name.

        :param str name: name of the CoordinateFrame.
        :return str vrt: the CoordinateFrame with the parsed name.
        """
        for vrt in self._vrts:
            if name == vrt.name:
                return vrt
        raise ValueError(name + " does not exist")

    def get_unapplied_vrt(self, name):
        """Returns the CoordinateFrame with the given name from the unapplied_vrt list.

        :param str name: name of the CoordinateFrame.
        :return str vrt: the CoordinateFrame with the parsed name.
        """
        for vrt in self._init_vrts:
            if name == vrt.name:
                return vrt
        raise ValueError(name + " does not exist")

    def get_downstream_connections(self, jump):
        """Finds the vrts that are related downstream from a jump.

        :param str jump: the jump name.
        :return: a list of vrt names.
        """
        vrts = []
        vrts.append(self._jumps[jump][1])
        dead_end = False
        while not dead_end:
            for counter, vrt_pair in enumerate(self._jumps.values(), 1):
                if vrt_pair[0] in vrts:
                    next_vrt = vrt_pair[1]
                    vrts.append(next_vrt)
                if counter == len(self._jumps.values()):
                    dead_end = True
        return vrts

    def get_master_jumps(self):
        """Gets the jumps that are master jumps.

        :return list: list of the names of the master jumps.
        """
        master_jumps = []
        for jumpgroup, jumps in self._jumpgroups.items():
            master_jumps.append(jumps[0])
        return master_jumps

    def is_dof_master_jumps(self):
        """Test if the jumps for which the dofs are set are also the master jumps (which it should be).

        :return bool: true/false.
        """
        for jump in self._dofs:
            if jump not in self.get_master_jumps():
                return False
        return True

    def make_symmetry_definition(self, use_stored_anchor=False, anchor_moved_resnums=0, headers=None, inherit_headers=True,
                                 use_current_vrt_positions=False):
        """Writes the symmetry definition to a python string."""
        symdef = ""
        # first write headers
        if inherit_headers:
            for k, v in self.headers.items():
                symdef += f"#{k}={v}\n"
        if headers is not None:
            for header in headers:
                symdef += f"#{header}" + "\n"
        symdef += "symmetry_name " + self.symmetry_name + "\n"
        symdef += "E = " + self.energies + "\n"
        symdef += "anchor_residue " + (str(self.actual_anchor_residue + anchor_moved_resnums) if use_stored_anchor else self.anchor) + "\n"
        if self.recenter:
            symdef += "recenter" + "\n"
        symdef += "virtual_coordinates_start" + "\n"
        # use either the potentially modified or the unapplied vrts to specifiy
        if use_current_vrt_positions or self._init_vrts is None:
            for vrt in self._vrts:
                symdef += str(vrt) + "\n"
        else:
            for vrt in self._init_vrts:
                symdef += str(vrt) + "\n"
        symdef += "virtual_coordinates_stop" + "\n"
        for name, connection in self._jumps.items():
            symdef += "connect_virtual " + name + " " + connection[0] + " " + connection[1] + "\n"
        for jump, dofs in self._dofs.items():
            symdef += "set_dof " + jump
            for dof in dofs:
                if dof[1] == "translation":
                    symdef += " " + dof[0]
                else:
                    symdef += " angle_" + dof[0]
                if dof[2] is not None:
                    symdef += "(" + str(dof[2]) + ")"
            symdef += "\n"
        for i, (name, jumps) in enumerate(self._jumpgroups.items(), 1):
            symdef += "set_jump_group " + name
            for jump in jumps:
                symdef += " " + jump
            if not i == len(self._jumpgroups):
                symdef += "\n"
        return symdef

    def output(self, name, headers=None, inherit_headers=True):
        """Prints the symmetry file to disk.

        :param str name: name given to file.
        """
        with open(name, 'w') as f:
            f.write(self.make_symmetry_definition(headers=headers, inherit_headers=inherit_headers))


    def get_dof_value(self, jumpname:str, dof:str, doftype:str):
        """Get the value for the dof"""
        for idof in self._dofs[jumpname]:
            if idof[0] == dof and idof[1] == doftype:
                return idof[2]
        raise ValueError(f"value for jumpname:{jumpname}, dof:{dof}, doftype:{doftype} not found!")

    def get_anchor_residue(self, pose):
        """Gets the anchor residue."""
        if self.anchor == "COM":
            return residue_center_of_mass(pose.conformation(), 1, pose.chain_end(1))
        else:
            raise NotImplementedError

    def read_from_pose_comment(self, pose, save_anchor_residue=False):
        """Reads a symmetry file from pose.

        :param pose: pose to read from. Must contain SYMMETRY in the datacache.
        :param save_anchor_residue: saves the actual anchor residue. Useful if the pose has changed after reading the symmetry for
        the first time, for example if actions needs to be taken on the asymmetric pose and the pose needs to be symmetrized again.
        :return: None
        """
        syminfo = pose.data().get_ptr(CacheableDataType.STRING_MAP).map()["SYMMETRY"].split(" | ")
        self.extract_symmetry_info(syminfo)
        self.reference_symmetric = self.is_reference_symmetry()
        if save_anchor_residue:
            self.actual_anchor_residue = self.get_anchor_residue(pose)
        self._set_init_vrts()

    def read_from_file(self, filename, check_for_reference_symmetry=True):
        """Reads a symmetry file from disk.

        :param filename: name of file to read from.
        :return: None
        """
        if isinstance(filename, StringIO):
            file = filename
        else:
            file = open(filename, 'r')
        self.extract_symmetry_info(file)
        if check_for_reference_symmetry:
            self.reference_symmetric = self.is_reference_symmetry()
        else:
            self.reference_symmetric = None
        self._set_init_vrts()

    @staticmethod
    def output_vrts_as_pdb(pose, filename):
        vrts = []
        for ri in range(1, pose.size() + 1):
            if pose.residue(ri).name() == "VRT":
                vrts.append(ri)
        string = "{:<6}{:>5}  {:<3}{:>4} {}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}        {:>4}{}\n"
        tmp_f = filename
        with open(tmp_f, "w") as f:
            i = 1
            for n, vri in enumerate(vrts, 1):
                o = np.array(pose.residue(vri).atom(1).xyz())
                x = np.array(pose.residue(vri).atom(2).xyz())
                y = np.array(pose.residue(vri).atom(3).xyz())
                # ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
                #                      ATOM     1  N    PRO    A      1  8.316 21.20 21.5  1.00 17.44      N
                f.write(string.format("HETATM", i, "O", "VRT", "X", n, o[0], o[1], o[2], 1.00, 0.00, "", ""))
                i +=1
                f.write(string.format("HETATM", i, "X", "VRT", "X", n, x[0], x[1], x[2], 1.00, 0.00, "", ""))
                i +=1
                f.write(string.format("HETATM", i, "Y", "VRT", "X", n, y[0], y[1], y[2], 1.00, 0.00, "", ""))
                i += 1
        # string = "{:<6}{:>5d}  {:<3}{:>4} {}{:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}        " \

        #          "{:>4}{}".format("HETATM", int(subunit.id), subunit.id, "gc", "P",
        #                           int(subunit.id), center_of_mass[0], center_of_mass[1],
        #                           center_of_mass[2], 1.00, 0.00,
        #                           "X", "")

    def extract_symmetry_info(self, info):
        for line in info:
            line = line.split()
            # headers
            if line[0][0] == "#":
                k, v = line[0][1:].split("=")
                self.headers[k] = v
            elif line[0] == "symmetry_name":
                self.symmetry_name = " ".join(line[1:])
            elif line[0] == "virtual_transforms_start":
                pass
                # using denovo script!!!
                # can we transform it into
            elif line[0] == "E":
                self.energies = " ".join(line[2:])
            elif line[0] == "anchor_residue":
                self.anchor = " ".join(line[1:])
            elif line[0] == "recenter":
                self.recenter = True
            elif line[0] == "xyz":
                vrt_name = line[1]
                x = np.array(line[2].split(","), dtype=np.float)
                y = np.array(line[3].split(","), dtype=np.float)
                z = np.cross(x, y)
                self.add_vrt((CoordinateFrame(vrt_name, x, y, z, np.array(line[4].split(","), dtype=np.float))))
            elif line[0] == "connect_virtual":
                self.add_jump(line[1], line[2], line[3])
            elif line[0] == "set_dof":
                for dof in line[2:]:
                    if dof[:5] == "angle":
                        degree = "rotation"
                        axes = dof[6]
                    else:
                        degree = "translation"
                        axes = dof[0]
                    if "(" in dof:
                        value = float(dof.split("(")[1].split(")")[0])
                    else:
                        value = None
                    self.add_dof(line[1], axes, degree, value)
            elif line[0] == "set_jump_group":
                self.add_jumpgroup(line[1], *line[2:])

    def visualize(self, apply_dofs=True, mark_jumps=True, ip="localhost", port="9123", suffix=""):
        """Visualizes the symmetry directly in PyMOL.

        :param bool apply_dofs: applies the translational and rotational degrees of freedom specified.
        :param bool mark_jumps: shows the jumps.
        :param str ip: the ip address of the machine where PyMOL is running.
        :param str ip: the port PyMOL is listening to.
        """
        cmd = xmlrpclib.ServerProxy(f'http://{ip}:{port}')
        cmd.do(self.__make_visualization_str(apply_dofs, mark_jumps, suffix=suffix))

    #    jump_name = jump_num_sym_dof(pose, jump_id)
    #         flexible_jump = pose.jump(jump_id)
    #         rot_euler = get_rotation_euler(flexible_jump)
    #         trans = get_translation(flexible_jump)
    #         for dof in range(1, 7):
    #             if symdof.allow_dof(dof):
    #                 if position_info.get(jump_name) == None:
    #                     position_info[jump_name] = {}
    #                 if dof == 1:
    #                     position_info[jump_name]["x"] = trans[0]
    #                 elif dof == 2:
    #                     position_info[jump_name]["y"] = trans[1]
    #                 elif dof == 3:
    #                     position_info[jump_name]["z"] = trans[2]
    #                 elif dof == 4:
    #                     position_info[jump_name]["angle_x"] = rot_euler[0]
    #                 elif dof == 5:
    #                     position_info[jump_name]["angle_y"] = rot_euler[1]
    #                 elif dof == 6:
    #                     position_info[jump_name]["angle_z"] = rot_euler[2]

    @staticmethod
    def get_jump(pose, jumpname):
        stub1 = Stub(pose.conformation().upstream_jump_stub(sym_dof_jump_num(pose, jumpname)))
        stub2 = Stub(pose.conformation().downstream_jump_stub(sym_dof_jump_num(pose, jumpname)))
        return Jump(stub1, stub2)

    def update_dofs_from_pose(self, pose, apply_dofs=False, assert_correct_overlap=True):
        """Updates the dofs from current dofs in the pose."""
        assert self.reference_symmetric, "Only works for reference symmetry."
        for jumpname, dofinfo in self._dofs.items():
            jump = self.get_jump(pose, jumpname)
            for dofname, doftype, _ in dofinfo:
                if doftype == "translation":
                    new_dofval = self.get_trans_from_jump(jump, dofname)
                else: # rotation
                    assert doftype == "rotation"
                    new_dofval = self.get_rot_from_jump(jump, "angle_" + dofname)
                for n, (t_dofname, t_doftype, _) in enumerate(self._dofs[jumpname]):
                    if t_dofname == dofname and t_doftype == doftype:
                        self._dofs[jumpname][n][2] = new_dofval
        if apply_dofs:
            self.apply_dofs()
        if assert_correct_overlap:
            self.vrts_overlap_with_pose(pose)

    def get_coordinateframes(self, apply_dofs=True):
        """Returns a list of the coordinates frames with or without applied dofs.

        :param apply_dofs: Applies the translational and rotational degrees of freedom specified in
                           the symmetry definition file.
        :return:
        """
        if apply_dofs and len(self._dofs) != 0:
            symmetry_setup = copy.deepcopy(self)
            symmetry_setup.apply_dofs()
        else:
            symmetry_setup = self
        return symmetry_setup._vrts

    def make_symmetric_pose(self, pose, use_stored_anchor=False, anchor_moved_resnums=0):
        """Symmetrizes the pose with the symmetrysetup (internal symmetry definition)."""
        s = SymmData()
        symdef = self.make_symmetry_definition(use_stored_anchor, anchor_moved_resnums)
        s.read_symmetry_data_from_stream(istringstream(symdef))
        setup = SetupForSymmetryMover(s)
        setup.apply(pose)

    def get_symmdata(self):
        """Constructs and returns a SymmData object based on the stored symmetry."""
        s = SymmData()
        s.read_symmetry_data_from_stream(istringstream(self.make_symmetry_definition()))
        return s

    @staticmethod
    def create_global_coordinateframe(name="VRTglobal"):
        """Returns a global CoordinateFrame object."""
        global_z = np.array([0, 0, 1])
        global_y = np.array([0, 1, 0])
        global_x = np.array([1, 0, 0])
        global_center = np.array([0, 0, 0])
        return CoordinateFrame(name, global_x, global_y, global_z, global_center)

    # def rotate_vrt_to_axis(self, vrtname, axis, angle):
    #     """Rotates the vrt with vrtname to the axis."""
    #     # rotation_matrix()
    #     pass

    def copy_vrt(self, vrtname, newname, move_origo=None, axis:str=None, mul=1, dir=-1):
        vrt_copy = copy.deepcopy(self.get_vrt(vrtname))
        if move_origo:
            if axis == "x":
                axis_v = vrt_copy._vrt_x
            elif axis == "y":
                axis_v = vrt_copy._vrt_y
            elif axis == "z":
                axis_v = vrt_copy._vrt_z
            else:
                raise ValueError(f"axis has to be either x, y or z not {axis}")
            vrt_copy._vrt_orig = self.add_along_vector(vrt_copy._vrt_orig, axis_v, mul, dir)
        vrt_copy.name = newname
        return vrt_copy

    def add_along_vector(self, v1, v2, mul=1, dir=-1):
        """move v1 along the direction of v2"""
        return v1 + (dir * v2 / np.linalg.norm(v2)) * mul

    @staticmethod
    def find_center_between_vtrs(setup, v1, v2):
        a = setup.get_vrt(v1).vrt_orig
        b = setup.get_vrt(v2).vrt_orig
        return (a + b) / 2

    def _set_init_vrts(self):
        """Everytime we apply the dofs we need to start from the initial vrts, therefore we save the vrts from the intial state."""
        self._init_vrts = copy.deepcopy(self._vrts)

    def _set_vrts_to_init_vrts(self):
        """Sets the current vrts to be of the initial vrts."""
        if self._init_vrts is None:
            raise ValueError("Initial vrts have not been set. These needs to be set before applying dofs. They can be set through: "
                             "1. Initializing the class with either a pose or a symdef file."
                             "2. Manual call to _set_init_vrts() if you are building a SymmetrySetup from scratch."
                             "3. Call to read_from_pose_comment() "
                             "4. Call to read_from_file()." )
        self._vrts = copy.deepcopy(self._init_vrts)

    def apply_dofs(self):
        """Applies the translational and rotational degrees of freedom specified in the symmetry definition file."""
        if not self.is_dof_master_jumps():
            raise ValueError("The DOFs that are set are not of the master jumps. Did you set the JUMPGROUP flag correctly?")
        # set the self.vrts to be the unapplied_vrts
        self._set_vrts_to_init_vrts()

        for jump_to_apply_dof_to, dofs in self._dofs.items():
            # find the jumpgroup the master jump degree of freedom belongs too
            for jumpgroup in self._jumpgroups.values():
                # check that the dof jump and the first entry (aka master jump) in the jumpgroup is the same
                # they could be out of order and therefore the first defined jump in set_dof might not be the
                # same in set_jump_group. Keep iterating until it is the same.
                # might need to check that none of them matches and output an error
                if jump_to_apply_dof_to == jumpgroup[0]:
                    break

            # apply dofs to all jumps in the jumpgroup
            for jump in jumpgroup:
                # find all downstream vrt names connected to the jump
                vrts_to_apply_dof_to = self.get_downstream_connections(jump)
                #  Find the reference vrt that the dofs should be applied from
                vrt_reference_name = self._jumps[jump][0]
                vrt_reference = self.get_vrt(vrt_reference_name)

                # now apply the dofs vrts_to_apply_dof_to
                for vrt_to_name in vrts_to_apply_dof_to:
                    if vrt_to_name == "SUBUNIT":
                        continue
                    vrt_to = self.get_vrt(vrt_to_name)
                    for dof in dofs:
                        if dof[2] is not None:
                            value = dof[2]
                        else:
                            continue
                        if dof[1] == "translation":
                            axis = dof[0]
                            if axis == 'x':
                                axis_to_apply_from = - vrt_reference.vrt_x  # minus because of Rosettas convention
                            elif axis == 'y':
                                axis_to_apply_from = vrt_reference.vrt_y
                            elif axis == 'z':
                                axis_to_apply_from = - vrt_reference.vrt_z  # minus because of Rosettas convention
                            vrt_to.vrt_orig = vrt_to.vrt_orig + axis_to_apply_from * value
                        elif dof[1] == "rotation":
                            axis = dof[0]
                            if axis == 'x':
                                axis_to_apply_from = - vrt_reference.vrt_x
                            elif axis == 'y':
                                axis_to_apply_from = vrt_reference.vrt_y
                            elif axis == 'z':
                                axis_to_apply_from = - vrt_reference.vrt_z
                            R = rotation_matrix(axis_to_apply_from, value)
                            vrt_to.vrt_x = rotate(vrt_to.vrt_x, R)
                            vrt_to.vrt_y = rotate(vrt_to.vrt_y, R)
                            vrt_to.vrt_z = rotate(vrt_to.vrt_z, R)
                            vrt_to._vrt_orig = rotate(vrt_to._vrt_orig - vrt_reference._vrt_orig, R) + vrt_reference._vrt_orig
                            # vrt_to.vrt_orig = rotate(vrt_to._vrt_orig, R)
                            # vrt_to.vrt_orig + axis_to_apply_from * value

    def __make_visualization_str(self, apply_dofs=True, mark_jumps=True, suffix=""):
        """Makes python script as a str that can either be printed to a file or use in PyMOL directly.

        :param bool apply_dofs: applies the translational and rotational degrees of freedom specified.
        :param bool mark_jumps: shows the jumps.
        """
        string = ""

        # CGO styles
        w = 0.06 * 3  # cylinder width
        w2 = 0.03  # cylinder width for jumps
        h = 0.06 * 3  # cone height
        d = 0.13 * 3  # cone base diameter
        r = 0.12 * 3  # sphere radius
        size = 0.01 * 6
        color = [1.0, 1.0, 1.0]  # color

        string += textwrap.dedent(
            '''
            # to produce vectors and spheres (VRT vectors and positions) in PDB
            from pymol.cgo import *
            from pymol import cmd
            from pymol.vfont import plain

            w = {0} # cylinder width
            w2 = {6} # cylinder width
            h = {1} # cone hight
            d = {2} # cone base diameter
            r = {3} # sphere radius
            size = {4}
            color= {5}

            # SPHERE, x, y, z,  radius 
            # CYLINDER, x1, y1, z1, x2, y2, z2, radius, red1, green1, blue1, red2, green2, blue2,
            # CONE, x0, y0, z0, x1, y1, z1, radius, r1, g1, b1, r1, g1, b1 

            ''').format(
            w, h, d, r, size, color, w2)

        # tmr you need to go through all the vrts that are related to a jump. not just the jumpgroup + make sure you
        # understand the how to change the objects in a list.

        # if user want to appy dofs and dofs are set then apply the dofs
        if apply_dofs and len(self._dofs) != 0:
            symmetry_setup = copy.deepcopy(self)
            symmetry_setup.apply_dofs()
        else:
            symmetry_setup = self

        # For text position
        p = 1.1  # arrow pointiness
        name_pos = 0.5
        x_pos = 1.2
        y_pos = 1.3
        z_pos = 1.2
        axes_norm = 0.3

        for counter, vrt in enumerate(symmetry_setup._vrts, 1):
            string += "obj{0} = [SPHERE, {1}, r," \
                      "CYLINDER, {1}, {2}, w, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0," \
                      "CYLINDER, {1}, {3}, w, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0," \
                      "CYLINDER, {1}, {4}, w, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, " \
                      "CONE, {2}, {5}, d, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0," \
                      "CONE, {3}, {6}, d, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0," \
                      "CONE, {4}, {7}, d, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,]".format(
                counter,
                ','.join(map(str, vrt._vrt_orig)),  # SPHERE
                ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * 3)),  # CYLINDER
                ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * 3)),  # CYLINDER
                ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * 3)),  # CYLINDER
                ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * p * 3)),  # CONE
                ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * p * 3)),  # CONE
                ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * p * 3)))  # CONE

            # TODO: after some refactoring some numbers are not important anymore
            string += textwrap.dedent(
                '''
                cyl_text(obj{0},plain,[{8}],'{12}',size, color, {13})
                cyl_text(obj{0},plain,[{9}],'x',size,color, {13})
                cyl_text(obj{0},plain,[{10}],'y',size,color, {13})
                cyl_text(obj{0},plain,[{11}],'z',size,color, {13})

                cmd.load_cgo(obj{0}, '{12}')
                '''.format(
                    counter,
                    ','.join(map(str, vrt._vrt_orig)),  # SPHERE
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * 3)),  # CYLINDER
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * 3)),  # CYLINDER
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * 3)),  # CYLINDER
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * p * 3)),  # CONE
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * p * 3)),  # CONE
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * p * 3)),  # CONE
                    ','.join(map(str, vrt._vrt_orig - vrt._vrt_y * name_pos * 1.5)),  # cyl_text
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * 3 * x_pos)),  # cyl_text
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * 3 * y_pos)),  # cyl_text
                    ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * 3 * z_pos)),  # cyl_text
                    vrt.name + suffix,  # name of vrt
                    'axes = ' + '[' + str((vrt._vrt_x * axes_norm).tolist()) + ',' + str(
                        (vrt._vrt_y * axes_norm).tolist()) + ',' + str((vrt._vrt_z * axes_norm).tolist()) + ']',
                    # cyl_text
                ))

        # mark jumps:
        if mark_jumps:

            # CGO styles
            w = 0.06  # cylinder width

            for counter2, (jump, vrts) in enumerate(symmetry_setup._jumps.items(), counter):
                vrt_from = symmetry_setup.get_vrt(vrts[0])
                if vrts[1] == "SUBUNIT":
                    continue
                vrt_to = symmetry_setup.get_vrt(vrts[1])
                string += textwrap.dedent((
                    """

                    obj{0} = [CYLINDER, {1},{2}, w2, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,]
                    cmd.load_cgo(obj{0}, '{3}')

                    """.format(
                        counter2,
                        ','.join(map(str, vrt_from._vrt_orig)),  # CYLINDER
                        ','.join(map(str, vrt_to._vrt_orig)),  # CYLINDER
                        jump,
                    )))

        return string

    def print_visualization(self, name, apply_dofs=True, mark_jumps=True):
        """Prints a python script that can be run in PyMOL to generate a visualization of the symmetric setup.

        :param str name: name given to the script
        :param bool apply_dofs: if to apply the degrees of freedom of the symmetric system or not.
        :param bool mark_jumps: shows the jumps.
        """

        file = open(name, 'w')
        file.write(self.__make_visualization_str(apply_dofs, mark_jumps))
        file.close()

    # def add_angle_vrt(self, pose):
    #     self.update_dofs_from_pose(pose)
    #     apose = self.make_asymmetric_pose(pose, dont_reset=["JUMPHFfold1111_subunit"])
    #     # get the attachmentpoint for
    #     vrt = CoordinateFrame("VRT_angle_z_controller"), [1, 0, 0], [0, 1, 0], [0, 0, 1], np.array(line[4].split(","), dtype=np.float)))
    #     self.add_vrt(vrt)
    #     ...
