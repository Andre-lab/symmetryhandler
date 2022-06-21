#!/usr/bin/env python3
# coding=utf-8
from .mathfunctions import rotation_matrix, rotation_matrix_from_vector_to_vector, rotate, vector_angle, vector_projection_on_subspace
import copy
import textwrap
import numpy as np
import xmlrpc.client as xmlrpclib
from io import StringIO
import warnings

# from shapedesign
from shapedesign.src.utilities.pose import get_position_info, dof_map
from pyrosetta.rosetta.std import istringstream
from pyrosetta.rosetta.core.conformation.symmetry import SymmData
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.pose.datacache import CacheableDataType
from shapedesign.src.utilities.pose import get_position_info, dof_map
from shapedesign.src.utilities.kinematics import set_jumpdof
from pyrosetta import Pose
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit
from pyrosetta.rosetta.core.conformation.symmetry import residue_center_of_mass

class UndefinedParameters(BaseException):
    """An exception that reports that the user has undefined parameters. Used in the CoordinateFrame class."""
    def __init__(self, message=None):
        self.message = message
        print(message)

class CoordinateFrame:
    """A coordinate system containing 4 virtual residues. It is used when modelling symmetry in Rosetta."""

    def __init__(self, name, x=None, y=None, z=None, orig=None, generate_z_from_xy=False):
        self.name = name

        if type(x) is list:
            self._vrt_x = np.array(x)
        else:
            self._vrt_x = x

        if type(y) is list:
            self._vrt_y = np.array(y)
        else:
            self._vrt_y = y

        if generate_z_from_xy:
            self.generate_z()
        else:
            if type(z) is list:
                self._vrt_z = np.array(z)
            else:
                self._vrt_z = z

        if type(orig) is list:
            self._vrt_orig = np.array(orig)
        else:
            self._vrt_orig = orig

    def generate_z(self):
        self.vrt_z = np.cross(self.vrt_x, self.vrt_y)

    def rotate(self, R):
        """Rotates the coordinate frame with the rotation matrix R"""
        self._vrt_orig = rotate(self._vrt_orig, R)
        self._vrt_x = rotate(self._vrt_x, R)
        self._vrt_y = rotate(self._vrt_y, R)
        self._vrt_z = rotate(self._vrt_z, R)

    @property
    def vrt_x(self):
        return self._vrt_x

    @vrt_x.setter
    def vrt_x(self, vrt_x):
        self._vrt_x = vrt_x

    @property
    def vrt_y(self):
        return self._vrt_y

    @vrt_y.setter
    def vrt_y(self, vrt_y):
        self._vrt_y = vrt_y

    @property
    def vrt_x(self):
        return self._vrt_x

    @vrt_x.setter
    def vrt_x(self, vrt_x):
        self._vrt_x = vrt_x

    @property
    def vrt_z(self):
        return self._vrt_z

    @vrt_z.setter
    def vrt_z(self, vrt_z):
        self._vrt_z = vrt_z

    @property
    def vrt_orig(self):
        return self._vrt_orig

    @vrt_orig.setter
    def vrt_orig(self, vrt_orig):
        self._vrt_orig = vrt_orig

    def __str__(self):
        str_rosetta_style = "xyz " + self.name
        for vrt in (self._vrt_x, self._vrt_y, self._vrt_orig):
            if vrt is None:  # do not print anything if these are not defined
                raise UndefinedParameters("parameters of the CoodinateFrame have not been defined")
            str_rosetta_style += " "
            for i, index in enumerate(vrt, 1):
                str_rosetta_style += "{:.6f}".format(index)
                if i != 3:
                    str_rosetta_style += ","
        return str_rosetta_style

class SymmetrySetup:
    """A symmetric setup that stores the symmetry information used internally in Rosetta."""

    def __init__(self, symmetry_name=None, anchor=None, energies=None, recenter=False, vrts=None, jumps=None,
                 jumpgroups=None, dofs=None):
        """Initializes the Symmetric setup.

        See this link to the Rosetta user page:
        https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry#symmetry-definitions

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
        self.anchor = anchor
        self.recenter = recenter
        if energies == None:
            self.energies = []
        else:
            self.energies = energies
        if jumps == None:
            self._jumps = {}
        else:
            self._jumps = jumps
        if jumpgroups == None:
            self._jumpgroups = {}
        else:
            self._jumpgroups = jumpgroups
        if dofs == None:
            self._dofs = {}
        else:
            self._dofs = dofs

        if vrts == None:
            self._vrts = []
        else:
            self._vrts = vrts

        #before with _ they become a part of the class and will be inherited every time an instance is creat
        # self._jumps = jumps
        # self._jumpgroups = jumpgroups
        # self._dofs = dofs
        # self._vrts = vrts

    def reset_all_dofs(self):
        """Resets all the values of the dofs."""
        for jn, vl in self._dofs.items():
            for v in vl:
                v[2] = 0.0

    def set_dof(self, jumpname, dof, doftype, value):
        """Sets the value of the dof.

        :param jumpname: Jumpname
        :param dof: dof (x,y,z)
        :param doftype: (translation, rotation)
        :param value: value to set for the dof
        :return:
        """
        for jn, vl in self._dofs.items():
            if jn == jumpname:
                for n, v in enumerate(vl):
                    t_dof, t_doftype, _ = v
                    if t_dof == dof and t_doftype == doftype:
                        self._dofs[jn][n][2] = value

    def add_vrt(self, vrt):
        """Adds a CoordinateFrame to this instance.

        :param CoordinateFrame vrt: name of the CoordinateFrame.
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

    def get_vrt_name(self, name):
        """Returns the CoordinateFrame with the given name.

        :param str name: name of the CoordinateFrame.
        :return str vrt: the CoordinateFrame with the parsed name.

        """
        for vrt in self._vrts:
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

    def make_symmetry_definition(self, use_stored_anchor=False, anchor_moved_resnums=0):
        """Writes the symmetry definition to a python string."""
        symdef = ""
        symdef += "symmetry_name " + self.symmetry_name + "\n"
        symdef += "E = " + self.energies + "\n"
        symdef += "anchor_residue " + (str(self.actual_anchor_residue + anchor_moved_resnums) if use_stored_anchor else self.anchor) + "\n"
        if self.recenter:
            symdef += "recenter" + "\n"
        symdef += "virtual_coordinates_start" + "\n"
        for vrt in self._vrts:
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

    def output(self, name):
        """Prints the symmetry file to disk.

        :param str name: name given to file.
        """
        with open(name, 'w') as f:
            f.write(self.make_symmetry_definition())

    def read_from_pose(self, pose, save_anchor_residue=False):
        """Reads a symmetry file from pose.

        :param pose: pose to read from. Must contain SYMMETRY in the datacache.
        :param save_anchor_residue: saves the actual anchor residue. Useful if the pose has changed after reading the symmetry for
        the first time, for example if actions needs to be taken on the asymmetric pose and the pose needs to be symmetrized again.
        :return: None
        """
        syminfo = pose.data().get_ptr(CacheableDataType.STRING_MAP).map()["SYMMETRY"].split(" | ")
        self.extract_symmetry_info(syminfo)
        if save_anchor_residue:
            self.actual_anchor_residue = residue_center_of_mass(pose.conformation(), 1, pose.chain_end(1))

    def read_from_file(self, name):
        """Reads a symmetry file from disk.

        :param name: name of file to read from.
        :return: None
        """
        if type(name) == str:
            file = open(name, 'r')
        else: # if io.StringIO
            file = name
        # extract
        # info = [l.replace("\n", "") for l in file.readlines()]
        self.extract_symmetry_info(file)

    def get_asymmetric_pose(self, pose, reset_dofs=True, dont_reset:list=None):
        if reset_dofs:
            pose = pose.clone() # will shadow the pose and pose will not change
            # set all degrees of freedom to 0
            for jump, dofs in get_position_info(pose, dictionary=True).items():
                if dont_reset and jump in dont_reset:
                    continue
                for dof, old_val in dofs.items():
                    set_jumpdof(pose, jump, dof_map[dof], 0)
        # create a new pose object and fill it with the asymmetric pose
        apose = Pose()
        extract_asymmetric_unit(pose, apose, False)
        return apose

    def extract_symmetry_info(self, info):
        for line in info:
            line = line.split()
            if line[0] == "symmetry_name":
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
                    value = dof.split("(")[1].split(")")[0]
                    self.add_dof(line[1], axes, degree, float(value))
            elif line[0] == "set_jump_group":
                self.add_jumpgroup(line[1], *line[2:])

    def visualize(self, apply_dofs=True, mark_jumps=True, ip="localhost", port="9123"):
        """Visualizes the symmetry directly in PyMOL.

        :param bool apply_dofs: applies the translational and rotational degrees of freedom specified.
        :param bool mark_jumps: shows the jumps.
        :param str ip: the ip address of the machine where PyMOL is running.
        :param str ip: the port PyMOL is listening to.
        """
        cmd = xmlrpclib.ServerProxy(f'http://{ip}:{port}')
        cmd.do(self.__make_visualization_str(apply_dofs, mark_jumps))

    def update_from_pose(self, pose):
        """Updates the dofs from current dofs in the pose."""
        # Get the set_dof lines in the symmetry file and see how the names have changed. Then modify the names
        position_info = get_position_info(pose, dictionary=True)
        # update the dofs
        for jumpname, dofinfo in position_info.items():
            dofs = self._dofs[jumpname]
            # fixme: change how this is accessed - this is a bit cluncky
            #  the self._dofs should have been a dict of dict
            new_dofs = []
            for dofname, dofval in dofinfo.items():
                 # condition means that the dofname contains 'angle' and is therefore a 'rotation'
                t = f"{'rotation' if len(dofname) > 1 else 'translation'}"
                d = f"{dofname.split('_')[1] if len(dofname) > 1 else dofname}"
                for dof in dofs:
                    if d == dof[0] and t == dof[1]:
                        new_dofs.append([*dof[:-1], dofval])
            self._dofs[jumpname] = new_dofs

    def get_coordinateframes(self, apply_dofs=True):
        """Returns a list of the coordinates frames with or without applied dofs.

        :param apply_dofs: Applies the translational and rotational degrees of freedom specified in
                           the symmetry definition file.
        :return:
        """
        if apply_dofs and len(self._dofs) != 0:
            symmetry_setup = copy.deepcopy(self)
            self.__apply_dofs(symmetry_setup)
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

    def rotations_to_2folds(self, pose, visualize=False, cmd=None):
        """Gets the rotation angles from the master subunit com to the two 2-fold symmetry axes.

        When looking at the icosahedral structure in PyMOL the neative angle_z rotation is to the left and the positive
        to the right.
        """
        symmetry_setup = copy.deepcopy(self)
        symmetry_setup.update_from_pose(pose)
        self.__apply_dofs(symmetry_setup)

        # main 5-fold vectors
        v_5fold_center_to_5fold_master_com = symmetry_setup.get_vrt_name("VRT5fold1111").vrt_orig - symmetry_setup.get_vrt_name("VRT5fold1").vrt_orig
        # v_5fold_center_to_5fold_slave2_com = symmetry_setup.get_vrt_name("VRT5fold1211").vrt_orig - symmetry_setup.get_vrt_name("VRT5fold1").vrt_orig
        # v_5fold_center_to_5fold_slave5_com = symmetry_setup.get_vrt_name("VRT5fold1511").vrt_orig - symmetry_setup.get_vrt_name("VRT5fold1").vrt_orig

        # other 5-fold vectors
        v_5fold_center_to_2fold_center = symmetry_setup.get_vrt_name("VRT2fold1").vrt_orig - symmetry_setup.get_vrt_name("VRT5fold1").vrt_orig
        # v_2fold_center_to_3fold_center = symmetry_setup.get_vrt_name("VRT3fold1").vrt_orig - symmetry_setup.get_vrt_name("VRT2fold1").vrt_orig
        v_5fold_center_to_3fold_center = symmetry_setup.get_vrt_name("VRT3fold1").vrt_orig - symmetry_setup.get_vrt_name("VRT5fold1").vrt_orig

        # project these onto subspace
        # NOTE: vector_projection_on_subspace assumes orthonormal vectors in subspace!!!
        #  Since the capsid is oriented in the z direction then x-y spane a plane slicing through it. We can use that.
        v_5fold_center_to_2fold_center_projected = vector_projection_on_subspace(v_5fold_center_to_2fold_center,
                                                                                 np.array([1, 0, 0]),
                                                                                 np.array([0, 1, 0]))
        v_5fold_center_to_3fold_center_projected = vector_projection_on_subspace(v_5fold_center_to_3fold_center,
                                                                                  np.array([1, 0, 0]),
                                                                                  np.array([0, 1, 0]))

        angle_to_nearest_2fold = vector_angle(v_5fold_center_to_5fold_master_com, v_5fold_center_to_2fold_center_projected)
        angle_to_furthest_2fold = vector_angle(v_5fold_center_to_5fold_master_com, v_5fold_center_to_3fold_center_projected)

        if visualize:
            cmd.do(f"pseudoatom v_5fold_center_to_5fold_master_com, pos={list(v_5fold_center_to_5fold_master_com)}")
            # cmd.do(f"v_5fold_center_to_5fold_master_com {symmetry_setup.get_vrt_name('VRT5fold1').vrt_orig}, {symmetry_setup.get_vrt_name("VRT5fold1111").vrt_orig})
            # cmd.do(f"pseudoatom v_5fold_center_to_5fold_slave2_com, pos={list(v_5fold_center_to_5fold_slave2_com)}")
            # cmd.do(f"pseudoatom v_5fold_center_to_5fold_slave5_com, pos={list(v_5fold_center_to_5fold_slave5_com)}")
            cmd.do(f"pseudoatom v_5fold_center_to_2fold_center, pos={list(v_5fold_center_to_2fold_center)}")
            cmd.do(f"pseudoatom v_5fold_center_to_3fold_center, pos={list(v_5fold_center_to_3fold_center)}")
            cmd.do(f"pseudoatom v_5fold_center_to_2fold_center_projected, pos={list(v_5fold_center_to_2fold_center_projected)}")
            cmd.do(f"pseudoatom v_5fold_center_to_3fold_center_projected, pos={list(v_5fold_center_to_3fold_center_projected)}")


        # todo: record this before hand
        # Now, the two 2 2-folds can be either right or left ot the master subunuit. To determine if they are right of left we can take the
        # cross product of one of them and sew how it aligns with the global z-axis. If it is -z (global axis), the 2-fold is to the left,
        # or the negative direction, while it is vica versa for the other. We arbitrarily pick the two-fold that is closest. This is
        # the one connected by VRT2fold1 to calculate the cross product from.
        z_value = np.cross(v_5fold_center_to_5fold_master_com, v_5fold_center_to_2fold_center_projected)[2]
        if z_value < 0: # the nearest twofold is to the left / in the negative angle_z rotation direction. [1stm case]
            max_negative_angle = -angle_to_nearest_2fold / 2
            max_positive_angle = angle_to_furthest_2fold / 2
        else: # the nearest twofold is to the right / in the positive angle_z rotation direction. [4v4m case]
            max_negative_angle = -angle_to_furthest_2fold / 2
            max_positive_angle = angle_to_nearest_2fold / 2

        total_angle = angle_to_nearest_2fold + angle_to_furthest_2fold
        if not np.isclose(total_angle, 72, rtol=1e-02):
            warnings.warn(f"The total 5-fold angle is not 72 degrees but {total_angle} degrees. Inaccuracies can arise "
                          f"if the crystal structure is not perfectly symmetrical (and hence the symmdef file might "
                          f"not be as well).")

        return max_negative_angle, max_positive_angle

    # TODO: If needed create a universal n-fold function
    def create_independent_4fold_symmetries(self, pose):
        """Creates independent symmetries for the 4-fold."""
        # chain 1-2 and chain 1-3
        symmetry_setup = copy.deepcopy(self)
        self.__apply_dofs(symmetry_setup)
        symmetry_setup.update_from_pose(pose)

        # what dofs are available in the old file

        chain1_2 = SymmetrySetup()
        chain1_2.read_from_file(
            StringIO(textwrap.dedent(f"""symmetry_name chain1_2 
            E = 4*VRT000111 + 4*(VRT000111:VRT000222)
            anchor_residue COM 
            recenter
            virtual_coordinates_start
            xyz VRTglobal 1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT0001 -1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT0002 -0.000000,-1.000000,0.000000 -1.000000,0.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT00011 -1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT00022 -0.000000,-1.000000,0.000000 -1.000000,0.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT000111 -1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT000222 -0.000000,-1.000000,0.000000 -1.000000,0.000000,0.000000 0.000000,0.000000,0.000000
            virtual_coordinates_stop
            connect_virtual JUMPG1 VRTglobal VRT0001
            connect_virtual JUMPG2 VRTglobal VRT0002
            connect_virtual JUMP1 VRT0001 VRT00011
            connect_virtual JUMP2 VRT0002 VRT00022
            connect_virtual JUMP11 VRT00011 VRT000111 
            connect_virtual JUMP22 VRT00022 VRT000222 
            connect_virtual JUMP111 VRT000111 SUBUNIT
            connect_virtual JUMP222 VRT000222 SUBUNIT
            set_dof JUMP1 x({symmetry_setup._dofs['JUMP1'][0][2]}) 
            set_dof JUMP11 angle_x({symmetry_setup._dofs['JUMP11'][0][2]}) angle_y({symmetry_setup._dofs['JUMP11'][1][2]}) angle_z({symmetry_setup._dofs['JUMP11'][2][2]})
            set_dof JUMP111 angle_x({symmetry_setup._dofs['JUMP111'][0][2]}) angle_y({symmetry_setup._dofs['JUMP111'][1][2]}) angle_z({symmetry_setup._dofs['JUMP111'][2][2]})
            set_jump_group MODIFIED_BASEJUMP1 JUMP1 JUMP2
            set_jump_group MODIFIED_BASEJUMP2 JUMP11 JUMP22
            set_jump_group MODIFIED_BASEJUMP3 JUMP111 JUMP222
            """)))

        chain1_3 = SymmetrySetup()
        chain1_3.read_from_file(
            StringIO(textwrap.dedent(f"""symmetry_name chain1_3 
            E = 4*VRT000111 + 2*(VRT000111:VRT000333)
            anchor_residue COM 
            recenter
            virtual_coordinates_start
            xyz VRTglobal 1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT0001 -1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT0003 1.000000,-0.000000,0.000000 -0.000000,-1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT00011 -1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT00033 1.000000,-0.000000,0.000000 -0.000000,-1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT000111 -1.000000,0.000000,0.000000 0.000000,1.000000,0.000000 0.000000,0.000000,0.000000
            xyz VRT000333 1.000000,-0.000000,0.000000 -0.000000,-1.000000,0.000000 0.000000,0.000000,0.000000
            virtual_coordinates_stop
            connect_virtual JUMPG1 VRTglobal VRT0001
            connect_virtual JUMPG3 VRTglobal VRT0003
            connect_virtual JUMP1 VRT0001 VRT00011
            connect_virtual JUMP3 VRT0003 VRT00033
            connect_virtual JUMP11 VRT00011 VRT000111 
            connect_virtual JUMP33 VRT00033 VRT000333 
            connect_virtual JUMP111 VRT000111 SUBUNIT
            connect_virtual JUMP333 VRT000333 SUBUNIT
            set_dof JUMP1 x({symmetry_setup._dofs['JUMP1'][0][2]}) 
            set_dof JUMP11 angle_x({symmetry_setup._dofs['JUMP11'][0][2]}) angle_y({symmetry_setup._dofs['JUMP11'][1][2]}) angle_z({symmetry_setup._dofs['JUMP11'][2][2]})
            set_dof JUMP111 angle_x({symmetry_setup._dofs['JUMP111'][0][2]}) angle_y({symmetry_setup._dofs['JUMP111'][1][2]}) angle_z({symmetry_setup._dofs['JUMP111'][2][2]})
            set_jump_group MODIFIED_BASEJUMP1 JUMP1 JUMP3
            set_jump_group MODIFIED_BASEJUMP2 JUMP11 JUMP33
            set_jump_group MODIFIED_BASEJUMP3 JUMP111 JUMP333
            """)))

        return chain1_2, chain1_3

    def create_independent_icosahedral_symmetries(self, pose):
        """Creates independent symmetries for the icosahedral 5-fold, 3-fold and two 2-folds."""
        symmetry_setup = copy.deepcopy(self)
        self.__apply_dofs(symmetry_setup)
        symmetry_setup.update_from_pose(pose)

        fold5 = SymmetrySetup()
        fold5.read_from_file(
        StringIO(textwrap.dedent(f"""symmetry_name 5fold
        E = 60*VRT5fold1111 + 60*(VRT5fold1111:VRT5fold1211) + 60*(VRT5fold1111:VRT5fold1311)
        anchor_residue COM
        virtual_coordinates_start
        {self.get_vrt_name("VRTglobal")}
        {self.get_vrt_name("VRT5fold")}
        {self.get_vrt_name("VRT5fold1")}
        {self.get_vrt_name("VRT5fold11")}
        {self.get_vrt_name("VRT5fold111")}
        {self.get_vrt_name("VRT5fold1111")}
        {self.get_vrt_name("VRT5fold12")}
        {self.get_vrt_name("VRT5fold121")}
        {self.get_vrt_name("VRT5fold1211")}
        {self.get_vrt_name("VRT5fold13")}
        {self.get_vrt_name("VRT5fold131")}
        {self.get_vrt_name("VRT5fold1311")}
        {self.get_vrt_name("VRT5fold14")}
        {self.get_vrt_name("VRT5fold141")}
        {self.get_vrt_name("VRT5fold1411")}
        {self.get_vrt_name("VRT5fold15")}
        {self.get_vrt_name("VRT5fold151")}
        {self.get_vrt_name("VRT5fold1511")}
        virtual_coordinates_stop  
        connect_virtual JUMP5fold VRTglobal VRT5fold
        connect_virtual JUMP5fold1 VRT5fold VRT5fold1
        connect_virtual JUMP5fold11 VRT5fold1 VRT5fold11
        connect_virtual JUMP5fold111 VRT5fold11 VRT5fold111
        connect_virtual JUMP5fold1111 VRT5fold111 VRT5fold1111
        connect_virtual JUMP5fold1111_subunit VRT5fold1111 SUBUNIT
        connect_virtual JUMP5fold12 VRT5fold1 VRT5fold12
        connect_virtual JUMP5fold121 VRT5fold12 VRT5fold121
        connect_virtual JUMP5fold1211 VRT5fold121 VRT5fold1211
        connect_virtual JUMP5fold1211_subunit VRT5fold1211 SUBUNIT
        connect_virtual JUMP5fold13 VRT5fold1 VRT5fold13
        connect_virtual JUMP5fold131 VRT5fold13 VRT5fold131
        connect_virtual JUMP5fold1311 VRT5fold131 VRT5fold1311
        connect_virtual JUMP5fold1311_subunit VRT5fold1311 SUBUNIT
        connect_virtual JUMP5fold14 VRT5fold1 VRT5fold14
        connect_virtual JUMP5fold141 VRT5fold14 VRT5fold141
        connect_virtual JUMP5fold1411 VRT5fold141 VRT5fold1411
        connect_virtual JUMP5fold1411_subunit VRT5fold1411 SUBUNIT
        connect_virtual JUMP5fold15 VRT5fold1 VRT5fold15
        connect_virtual JUMP5fold151 VRT5fold15 VRT5fold151
        connect_virtual JUMP5fold1511 VRT5fold151 VRT5fold1511
        connect_virtual JUMP5fold1511_subunit VRT5fold1511 SUBUNIT
        set_dof JUMP5fold1 z({symmetry_setup._dofs['JUMP5fold1'][0][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1'][1][2]})
        set_dof JUMP5fold111 x({symmetry_setup._dofs['JUMP5fold111'][0][2]})
        set_dof JUMP5fold1111 angle_x({symmetry_setup._dofs['JUMP5fold1111'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111'][2][2]})
        set_dof JUMP5fold1111_subunit angle_x({symmetry_setup._dofs['JUMP5fold1111_subunit'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111_subunit'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111_subunit'][2][2]})
        set_jump_group JUMPGROUP1 JUMP5fold1 
        set_jump_group JUMPGROUP2 JUMP5fold111 JUMP5fold121 JUMP5fold131 JUMP5fold141 JUMP5fold151 
        set_jump_group JUMPGROUP3 JUMP5fold1111 JUMP5fold1211 JUMP5fold1311 JUMP5fold1411 JUMP5fold1511 
        set_jump_group JUMPGROUP4 JUMP5fold1111_subunit JUMP5fold1211_subunit JUMP5fold1311_subunit JUMP5fold1411_subunit JUMP5fold1511_subunit 
        """)))

        # TODO: change the symmetry so that depending on if it is 4v4m or 1stm different symmetries have to be used

        fold3 = SymmetrySetup()
        fold3.read_from_file(
        StringIO(textwrap.dedent(f"""symmetry_name 3fold
        E = 60*VRT5fold1111 + 60*(VRT5fold1111:VRT3fold1111)
        anchor_residue COM
        virtual_coordinates_start
        {self.get_vrt_name("VRTglobal")}
        {self.get_vrt_name("VRT5fold")}
        {self.get_vrt_name("VRT5fold1")}
        {self.get_vrt_name("VRT5fold11")}
        {self.get_vrt_name("VRT5fold111")}
        {self.get_vrt_name("VRT5fold1111")}
        {self.get_vrt_name("VRT3fold")}
        {self.get_vrt_name("VRT3fold1")}
        {self.get_vrt_name("VRT3fold11")}
        {self.get_vrt_name("VRT3fold111")}
        {self.get_vrt_name("VRT3fold1111")}
        {self.get_vrt_name("VRT2fold")}
        {self.get_vrt_name("VRT2fold1")}
        {self.get_vrt_name("VRT2fold12")}
        {self.get_vrt_name("VRT2fold121")}
        {self.get_vrt_name("VRT2fold1211")}
        virtual_coordinates_stop
        connect_virtual JUMP5fold VRTglobal VRT5fold
        connect_virtual JUMP5fold1 VRT5fold VRT5fold1
        connect_virtual JUMP5fold11 VRT5fold1 VRT5fold11
        connect_virtual JUMP5fold111 VRT5fold11 VRT5fold111
        connect_virtual JUMP5fold1111 VRT5fold111 VRT5fold1111
        connect_virtual JUMP5fold1111_subunit VRT5fold1111 SUBUNIT
        connect_virtual JUMP3fold VRTglobal VRT3fold
        connect_virtual JUMP3fold1 VRT3fold VRT3fold1
        connect_virtual JUMP3fold11 VRT3fold1 VRT3fold11
        connect_virtual JUMP3fold111 VRT3fold11 VRT3fold111
        connect_virtual JUMP3fold1111 VRT3fold111 VRT3fold1111
        connect_virtual JUMP3fold1111_subunit VRT3fold1111 SUBUNIT
        connect_virtual JUMP2fold VRTglobal VRT2fold
        connect_virtual JUMP2fold1 VRT2fold VRT2fold1
        connect_virtual JUMP2fold12 VRT2fold1 VRT2fold12
        connect_virtual JUMP2fold121 VRT2fold12 VRT2fold121
        connect_virtual JUMP2fold1211 VRT2fold121 VRT2fold1211
        connect_virtual JUMP2fold1211_subunit VRT2fold1211 SUBUNIT
        set_dof JUMP5fold1 z({symmetry_setup._dofs['JUMP5fold1'][0][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1'][1][2]})
        set_dof JUMP5fold111 x({symmetry_setup._dofs['JUMP5fold111'][0][2]})
        set_dof JUMP5fold1111 angle_x({symmetry_setup._dofs['JUMP5fold1111'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111'][2][2]})
        set_dof JUMP5fold1111_subunit angle_x({symmetry_setup._dofs['JUMP5fold1111_subunit'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111_subunit'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111_subunit'][2][2]})
        set_jump_group JUMPGROUP1 JUMP5fold1 JUMP3fold1 JUMP2fold1
        set_jump_group JUMPGROUP2 JUMP5fold111 JUMP3fold111  JUMP2fold121
        set_jump_group JUMPGROUP3 JUMP5fold1111 JUMP3fold1111  JUMP2fold1211
        set_jump_group JUMPGROUP4 JUMP5fold1111_subunit JUMP3fold1111_subunit JUMP2fold1211_subunit
        """)))

        fold2_1 = SymmetrySetup()
        fold2_1.read_from_file(
        StringIO(textwrap.dedent(f"""symmetry_name 2fold_1
        E = 60*VRT5fold1111 + 30*(VRT5fold1111:VRT2fold1111)
        anchor_residue COM
        virtual_coordinates_start
        {self.get_vrt_name("VRTglobal")}
        {self.get_vrt_name("VRT5fold")}
        {self.get_vrt_name("VRT5fold1")}
        {self.get_vrt_name("VRT5fold11")}
        {self.get_vrt_name("VRT5fold111")}
        {self.get_vrt_name("VRT5fold1111")}
        {self.get_vrt_name("VRT2fold")}
        {self.get_vrt_name("VRT2fold1")}
        {self.get_vrt_name("VRT2fold11")}
        {self.get_vrt_name("VRT2fold111")}
        {self.get_vrt_name("VRT2fold1111")}
        virtual_coordinates_stop
        connect_virtual JUMP5fold VRTglobal VRT5fold
        connect_virtual JUMP5fold1 VRT5fold VRT5fold1
        connect_virtual JUMP5fold11 VRT5fold1 VRT5fold11
        connect_virtual JUMP5fold111 VRT5fold11 VRT5fold111
        connect_virtual JUMP5fold1111 VRT5fold111 VRT5fold1111
        connect_virtual JUMP5fold1111_subunit VRT5fold1111 SUBUNIT
        connect_virtual JUMP2fold VRTglobal VRT2fold
        connect_virtual JUMP2fold1 VRT2fold VRT2fold1
        connect_virtual JUMP2fold11 VRT2fold1 VRT2fold11
        connect_virtual JUMP2fold111 VRT2fold11 VRT2fold111
        connect_virtual JUMP2fold1111 VRT2fold111 VRT2fold1111
        connect_virtual JUMP2fold1111_subunit VRT2fold1111 SUBUNIT
        set_dof JUMP5fold1 z({symmetry_setup._dofs['JUMP5fold1'][0][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1'][1][2]})
        set_dof JUMP5fold111 x({symmetry_setup._dofs['JUMP5fold111'][0][2]})
        set_dof JUMP5fold1111 angle_x({symmetry_setup._dofs['JUMP5fold1111'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111'][2][2]})
        set_dof JUMP5fold1111_subunit angle_x({symmetry_setup._dofs['JUMP5fold1111_subunit'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111_subunit'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111_subunit'][2][2]})
        set_jump_group JUMPGROUP1 JUMP5fold1 JUMP2fold1
        set_jump_group JUMPGROUP2 JUMP5fold111 JUMP2fold111 
        set_jump_group JUMPGROUP3 JUMP5fold1111 JUMP2fold1111 
        set_jump_group JUMPGROUP4 JUMP5fold1111_subunit JUMP2fold1111_subunit 
        """)))

        fold2_2 = SymmetrySetup()
        fold2_2.read_from_file(
        StringIO(textwrap.dedent(f"""symmetry_name fold2_2 
        E = 60*VRT5fold1111 + 30*(VRT5fold1111:VRT3fold1211)
        anchor_residue COM
        virtual_coordinates_start
        {self.get_vrt_name("VRTglobal")}
        {self.get_vrt_name("VRT5fold")}
        {self.get_vrt_name("VRT5fold1")}
        {self.get_vrt_name("VRT5fold11")}
        {self.get_vrt_name("VRT5fold111")}
        {self.get_vrt_name("VRT5fold1111")}
        {self.get_vrt_name("VRT3fold")}
        {self.get_vrt_name("VRT3fold1")}
        {self.get_vrt_name("VRT3fold12")}
        {self.get_vrt_name("VRT3fold121")}
        {self.get_vrt_name("VRT3fold1211")}
        virtual_coordinates_stop
        connect_virtual JUMP5fold VRTglobal VRT5fold
        connect_virtual JUMP5fold1 VRT5fold VRT5fold1
        connect_virtual JUMP5fold11 VRT5fold1 VRT5fold11
        connect_virtual JUMP5fold111 VRT5fold11 VRT5fold111
        connect_virtual JUMP5fold1111 VRT5fold111 VRT5fold1111
        connect_virtual JUMP5fold1111_subunit VRT5fold1111 SUBUNIT
        connect_virtual JUMP3fold VRTglobal VRT3fold
        connect_virtual JUMP3fold1 VRT3fold VRT3fold1
        connect_virtual JUMP3fold12 VRT3fold1 VRT3fold12
        connect_virtual JUMP3fold121 VRT3fold12 VRT3fold121
        connect_virtual JUMP3fold1211 VRT3fold121 VRT3fold1211
        connect_virtual JUMP3fold1211_subunit VRT3fold1211 SUBUNIT
        set_dof JUMP5fold1 z({symmetry_setup._dofs['JUMP5fold1'][0][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1'][1][2]})
        set_dof JUMP5fold111 x({symmetry_setup._dofs['JUMP5fold111'][0][2]})
        set_dof JUMP5fold1111 angle_x({symmetry_setup._dofs['JUMP5fold1111'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111'][2][2]})
        set_dof JUMP5fold1111_subunit angle_x({symmetry_setup._dofs['JUMP5fold1111_subunit'][0][2]}) angle_y({symmetry_setup._dofs['JUMP5fold1111_subunit'][1][2]}) angle_z({symmetry_setup._dofs['JUMP5fold1111_subunit'][2][2]})
        set_jump_group JUMPGROUP1 JUMP5fold1 JUMP3fold1
        set_jump_group JUMPGROUP2 JUMP5fold111 JUMP3fold121
        set_jump_group JUMPGROUP3 JUMP5fold1111 JUMP3fold1211
        set_jump_group JUMPGROUP4 JUMP5fold1111_subunit JUMP3fold1211_subunit
        """)))

        # setup_3fold = SymmetrySetup("3fold")
        # vrtglobal = symmetry_setup.get_vrt_name("VRTglobal")
        # center_of_3fold = np.array([symmetry_setup.get_vrt_name(vrt).vrt_orig for vrt in ("VRT5fold1111", "VRT3fold1111", "VRT2fold1111")]).sum(axis=1) / 3
        # rotation_to_3fold = rotation_matrix_from_vector_to_vector(vrtglobal.vrt_orig, center_of_3fold)
        # vrt3fold = copy.deepcopy(vrtglobal).rotate(rotation_to_3fold)

        return fold5, fold3, fold2_1, fold2_2

    def get_icosahedral_boundary_box(self):
        pass # get how much you can rotate/translate

    def __apply_dofs(self, symmetry_setup):
        """Applies the translational and rotational degrees of freedom specified in the symmetry definition file."""
        # checks that the degrees of freedom set are from the master jumps
        if not symmetry_setup.is_dof_master_jumps():
            raise ValueError(textwrap.dedent("""
                                    The degrees of freedom set are not of the master jumps.
                                    The results will not make sense if not.
                                    """))
        for jump_to_apply_dof_to, dofs in symmetry_setup._dofs.items():
            # find the jumpgroup the master jump degree of freedom belongs too
            for jumpgroup in symmetry_setup._jumpgroups.values():
                # check that the dof jump and the first entry (aka master jump) in the jumpgroup is the same
                # they could be out of order and therefore the first defined jump in set_dof might not be the
                # same in set_jump_group. Keep iterating until it is the same.
                # might need to check that none of them matches and output an error
                if jump_to_apply_dof_to == jumpgroup[0]:
                    break

            # apply dofs to all jumps in the jumpgroup
            for jump in jumpgroup:
                # find all downstream vrt names connected to the jump
                vrts_to_apply_dof_to = symmetry_setup.get_downstream_connections(jump)
                #  Find the reference vrt that the dofs should be applied from
                vrt_reference_name = symmetry_setup._jumps[jump][0]
                vrt_reference = symmetry_setup.get_vrt_name(vrt_reference_name)

                # now apply the dofs vrts_to_apply_dof_to
                for vrt_to_name in vrts_to_apply_dof_to:
                    if vrt_to_name == "SUBUNIT":
                        continue
                    vrt_to = symmetry_setup.get_vrt_name(vrt_to_name)
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
                                axis_to_apply_from = - vrt_reference.vrt_x  # minus because of Rosettas convention (i think it applies to rot)
                            elif axis == 'y':
                                axis_to_apply_from = vrt_reference.vrt_y
                            elif axis == 'z':
                                axis_to_apply_from = vrt_reference.vrt_z  # minus because of Rosettas convention(i think it applies to rot)
                            R = rotation_matrix(axis_to_apply_from, value)
                            vrt_to.vrt_x = rotate(vrt_to.vrt_x, R)
                            vrt_to.vrt_y = rotate(vrt_to.vrt_y, R)
                            vrt_to.vrt_z = rotate(vrt_to.vrt_z, R)

    def __make_visualization_str(self, apply_dofs=True, mark_jumps=True):
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
            symmetry_setup = copy.deepcopy(
                self)  # is this important? yes beacuse we are applying dofs now to the symmetry_setup
            # checks that the degrees of freedom set are from the master jumps
            self.__apply_dofs(symmetry_setup)
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
                    vrt.name,  # name of vrt
                    'axes = ' + '[' + str((vrt._vrt_x * axes_norm).tolist()) + ',' + str(
                        (vrt._vrt_y * axes_norm).tolist()) + ',' + str((vrt._vrt_z * axes_norm).tolist()) + ']',
                    # cyl_text
                ))

        # mark jumps:
        if mark_jumps:

            # CGO styles
            w = 0.06  # cylinder width

            for counter2, (jump, vrts) in enumerate(symmetry_setup._jumps.items(), counter):
                vrt_from = symmetry_setup.get_vrt_name(vrts[0])
                if vrts[1] == "SUBUNIT":
                    continue
                vrt_to = symmetry_setup.get_vrt_name(vrts[1])
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
