#!/usr/bin/env python3
# coding=utf-8
from lib.mathfunctions import rotation_matrix
from lib.mathfunctions import rotate
import copy
import textwrap
import numpy as np
import xmlrpc.client as xmlrpclib

class UndefinedParameters(BaseException):
    """An exception that reports that the user has undefined parameters. Used in the CoordinateFrame class."""
    def __init__(self, message=None):
        self.message = message
        print(message)

class CoordinateFrame:
    """A coordinate system containing 4 virtual residues. It is used when modelling symmetry in Rosetta."""

    def __init__(self, name, x=None, y=None, z=None, orig=None):
        self.name = name

        if type(x) is list:
            self._vrt_x = np.ndarray(x)
        else:
            self._vrt_x = x

        if type(y) is list:
            self._vrt_y = np.ndarray(y)
        else:
            self._vrt_y = y

        if type(z) is list:
            self._vrt_z = np.ndarray(z)
        else:
            self._vrt_z = z

        if type(orig) is list:
            self._vrt_orig = np.ndarray(orig)
        else:
            self._vrt_orig = orig

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

    def rosetta_repr(self):
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

    def output(self, name):
        """Prints the symmetry file to disk.

        :param str name: name given to file.

        """
        file = open(name, 'w')
        file.write("symmetry_name " + self.symmetry_name + "\n")
        file.write("E = " + self.energies + "\n")
        file.write("anchor_residue " + self.anchor + "\n")
        if self.recenter:
            file.write("recenter" + "\n")
        file.write("virtual_coordinates_start" + "\n")
        for vrt in self._vrts:
            file.write(vrt.rosetta_repr() + "\n")
        file.write("virtual_coordinates_stop" + "\n")
        for name, connection in self._jumps.items():
            file.write("connect_virtual " + name + " " + connection[0] + " " + connection[1] + "\n")
        for jump, dofs in self._dofs.items():
            file.write("set_dof " + jump)
            for dof in dofs:
                if dof[1] == "translation":
                    file.write(" " + dof[0])
                else:
                    file.write(" angle_" + dof[0])
                if dof[2] is not None:
                    file.write("(" + str(dof[2]) + ")")
            file.write("\n")
        for i, (name, jumps) in enumerate(self._jumpgroups.items(), 1):
            file.write("set_jump_group " + name)
            for jump in jumps:
                file.write(" " + jump)
            if not i == len(self._jumpgroups):
                file.write("\n")
        file.close()

    def read_from_file(self, name):
        """Reads a symmetry file from disk.

        :param name: name of file to read from.
        :return: None
        """
        file = open(name, 'r')
        for line in file:
            line = line.split()
            if line[0] == "symmetry_name":
                self.symmetry_name = " ".join(line[1:])
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
                z = np.cross(x,y)
                self.add_vrt((CoordinateFrame(vrt_name, x, y, z, np.array(line[4].split(","),dtype=np.float ))))
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
                                    axis_to_apply_from = - vrt_reference.vrt_z  # minus because of Rosettas convention(i think it applies to rot)
                                R = rotation_matrix(axis_to_apply_from, value)
                                vrt_to.vrt_x = rotate(vrt_to.vrt_x, R)
                                vrt_to.vrt_y = rotate(vrt_to.vrt_y, R)
                                vrt_to.vrt_z = rotate(vrt_to.vrt_z, R)
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

        # # CGO styles
        # w = 0.06 * 3 # cylinder width
        # w2 = 0.03 # cylinder width for jumps
        # h = 0.06 * 3 # cone height
        # d = 0.13 * 3 # cone base diameter
        # r = 0.12 * 3 # sphere radius
        # size = 0.01 * 6
        # color = [1.0, 1.0, 1.0]  # color
        #
        # file.write(textwrap.dedent(
        #     '''
        #     # to produce vectors and spheres (VRT vectors and positions) in PDB
        #     from pymol.cgo import *
        #     from pymol import cmd
        #     from pymol.vfont import plain
        #
        #     w = {0} # cylinder width
        #     w2 = {6} # cylinder width
        #     h = {1} # cone hight
        #     d = {2} # cone base diameter
        #     r = {3} # sphere radius
        #     size = {4}
        #     color= {5}
        #
        #     # SPHERE, x, y, z,  radius
        #     # CYLINDER, x1, y1, z1, x2, y2, z2, radius, red1, green1, blue1, red2, green2, blue2,
        #     # CONE, x0, y0, z0, x1, y1, z1, radius, r1, g1, b1, r1, g1, b1
        #
        #     ''').format(
        #     w, h, d, r, size, color, w2))
        #
        # # tmr you need to go through all the vrts that are related to a jump. not just the jumpgroup + make sure you
        # # understand the how to change the objects in a list.
        #
        # # if user want to appy dofs and dofs are set then apply the dofs
        # if apply_dofs and len(self._dofs) != 0:
        #     symmetry_setup = copy.deepcopy(self) # is this important? yes beacuse we are applying dofs now to the symmetry_setup
        #     # checks that the degrees of freedom set are from the master jumps
        #     if not symmetry_setup.is_dof_master_jumps():
        #         raise ValueError(textwrap.dedent("""
        #                          The degrees of freedom set are not of the master jumps.
        #                          The results will not make sense if not.
        #                          """))
        #     for jump_to_apply_dof_to, dofs in symmetry_setup._dofs.items():
        #         # find the jumpgroup the master jump degree of freedom belongs too
        #         for jumpgroup in symmetry_setup._jumpgroups.values():
        #             # check that the dof jump and the first entry (aka master jump) in the jumpgroup is the same
        #             # they could be out of order and therefore the first defined jump in set_dof might not be the
        #             # same in set_jump_group. Keep iterating until it is the same.
        #             # might need to check that none of them matches and output an error
        #             if jump_to_apply_dof_to == jumpgroup[0]:
        #                 break
        #
        #         # apply dofs to all jumps in the jumpgroup
        #         for jump in jumpgroup:
        #             # find all downstream vrt names connected to the jump
        #             vrts_to_apply_dof_to = symmetry_setup.get_downstream_connections(jump)
        #             #  Find the reference vrt that the dofs should be applied from
        #             vrt_reference_name = symmetry_setup._jumps[jump][0]
        #             vrt_reference = symmetry_setup.get_vrt_name(vrt_reference_name)
        #
        #             # now apply the dofs vrts_to_apply_dof_to
        #             for vrt_to_name in vrts_to_apply_dof_to:
        #                 if vrt_to_name == "SUBUNIT":
        #                     continue
        #                 vrt_to = symmetry_setup.get_vrt_name(vrt_to_name)
        #                 for dof in dofs:
        #                     if dof[2] is not None:
        #                         value = dof[2]
        #                     else:
        #                         continue
        #                     if dof[1] == "translation":
        #                         axis = dof[0]
        #                         if axis == 'x':
        #                             axis_to_apply_from = - vrt_reference.vrt_x # minus because of Rosettas convention
        #                         elif axis == 'y':
        #                             axis_to_apply_from = vrt_reference.vrt_y
        #                         elif axis == 'z':
        #                             axis_to_apply_from = - vrt_reference.vrt_z # minus because of Rosettas convention
        #                         vrt_to.vrt_orig = vrt_to.vrt_orig + axis_to_apply_from * value
        #                     elif dof[1] == "rotation":
        #                         axis = dof[0]
        #                         if axis == 'x':
        #                             axis_to_apply_from = - vrt_reference.vrt_x # minus because of Rosettas convention (i think it applies to rot)
        #                         elif axis == 'y':
        #                             axis_to_apply_from = vrt_reference.vrt_y
        #                         elif axis == 'z':
        #                             axis_to_apply_from = - vrt_reference.vrt_z # minus because of Rosettas convention(i think it applies to rot)
        #                         R = rotation_matrix(axis_to_apply_from, value)
        #                         vrt_to.vrt_x = rotate(vrt_to.vrt_x, R)
        #                         vrt_to.vrt_y = rotate(vrt_to.vrt_y, R)
        #                         vrt_to.vrt_z = rotate(vrt_to.vrt_z, R)
        # else:
        #     symmetry_setup = self
        #
        # # For text position
        # p = 1.1  # arrow pointiness
        # name_pos = 0.5
        # x_pos = 1.2
        # y_pos = 1.3
        # z_pos = 1.2
        # axes_norm = 0.3
        #
        # for counter, vrt in enumerate(symmetry_setup._vrts, 1):
        #     file.write(textwrap.dedent(
        #         '''
        #
        #         obj{0} = [
        #         SPHERE, {1}, r,
        #         CYLINDER, {1},{2}, w, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        #         CYLINDER, {1},{3}, w, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        #         CYLINDER, {1},{4}, w, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        #         CONE, {2},{5}, d, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        #         CONE, {3},{6}, d, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        #         CONE, {4},{7}, d, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        #            ]
        #
        #         cyl_text(obj{0},plain,[{8}],'{12}',size, color, {13})
        #         cyl_text(obj{0},plain,[{9}],'x',size,color, {13})
        #         cyl_text(obj{0},plain,[{10}],'y',size,color, {13})
        #         cyl_text(obj{0},plain,[{11}],'z',size,color, {13})
        #
        #         cmd.load_cgo(obj{0}, '{12}')
        #
        #         '''.format(
        #             counter,
        #             ','.join(map(str, vrt._vrt_orig)),  # SPHERE
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * 3)),  # CYLINDER
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * 3)),  # CYLINDER
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * 3)),  # CYLINDER
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * p * 3)),  # CONE
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * p * 3)),  # CONE
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * p * 3)),  # CONE
        #             ','.join(map(str, vrt._vrt_orig - vrt._vrt_y * name_pos * 1.5)),  # cyl_text
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_x * 3 * x_pos)),  # cyl_text
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_y * 3 * y_pos)),  # cyl_text
        #             ','.join(map(str, vrt._vrt_orig + vrt._vrt_z * 3 * z_pos)),  # cyl_text
        #             vrt.name,  # name of vrt
        #             'axes = ' + '[' + str((vrt._vrt_x * axes_norm).tolist()) + ',' + str(
        #                 (vrt._vrt_y * axes_norm).tolist()) + ',' + str((vrt._vrt_z * axes_norm).tolist()) + ']',
        #             # cyl_text
        #     )))
        #
        # #mark jumps:
        # if mark_jumps:
        #
        #     # CGO styles
        #     w = 0.06  # cylinder width
        #
        #     for counter2, (jump, vrts) in enumerate(symmetry_setup._jumps.items(), counter):
        #         vrt_from = symmetry_setup.get_vrt_name(vrts[0])
        #         if vrts[1] == "SUBUNIT":
        #             continue
        #         vrt_to = symmetry_setup.get_vrt_name(vrts[1])
        #         file.write(textwrap.dedent((
        #             """
        #
        #             obj{0} = [CYLINDER, {1},{2}, w2, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,]
        #             cmd.load_cgo(obj{0}, '{3}')
        #
        #             """.format(
        #                 counter2,
        #                 ','.join(map(str, vrt_from._vrt_orig)),  # CYLINDER
        #                 ','.join(map(str, vrt_to._vrt_orig)),  # CYLINDER
        #                 jump,
        #     ))))

        file.close()
