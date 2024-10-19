import shutil
import subprocess
import math
import tempfile
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from pyrosetta import pose_from_file
from pyrosetta.rosetta.core.conformation.symmetry import residue_center_of_mass
import symmetryhandler
from symmetryhandler.coordinateframe import CoordinateFrame
from symmetryhandler.mathfunctions import rotation_matrix_from_vector_to_vector as rvv, vector_projection, normalize, \
    vector_angle, rotation_matrix
from symmetryhandler.symmetrysetup import SymmetrySetup
from Bio.PDB import MMCIFParser, PDBIO, Select
import gzip

def cif2pdb(ciffile, structid, outname):
    parser = MMCIFParser()
    structure = parser.get_structure(structid,  gzip.open(ciffile,"rt"))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(outname))

class ChainSelect(Select):
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        if chain.get_id() in self.chain_letters:
            return 1
        else:
            return 0

class SymmetryMaker:

    make_symmdef_file = Path(symmetryhandler.__file__).parent.joinpath("scripts/make_symmdef_file.pl")
    print(make_symmdef_file)
    # assert make_symmdef_file.exists()

    def __init__(self, struct, point_group, chains, out, a, i=None, b=None, align_struct=None):
        self.struct_file = struct
        self.point_group = point_group
        self.out = out
        self.a = " ".join(a)
        self.i = " ".join(i)
        self.b = b
        self.chains_type = chains
        if align_struct is not None:
            parser = PDBParser(QUIET=True)
            self.align_struct = parser.get_structure("align_struct", align_struct.split("-")[0])
            self.align_chain = align_struct.split("-")[1]
        else:
            self.align_struct, self.align_chain  = None, None

    def get_center_of_mass_of_chain(self, chain):
        # Initialize total mass and center of mass vector
        total_mass = 0.0
        center_of_mass = np.zeros(3)

        # Iterate over atoms in the chain
        for atom in chain.get_atoms():
            mass = atom.mass  # Get the mass of the atom
            position = atom.coord  # Get the coordinates of the atom
            center_of_mass += mass * position
            total_mass += mass

        # Normalize by total mass to get the center of mass
        center_of_mass /= total_mass
        return center_of_mass

    def extract_struct_atoms(self, struct, atom_selection="CA"):
        """Extract atom coordinates based on atom_selection, e.g., C-alpha (CA)."""
        atoms = []
        for model in struct:
            for chain in model:
                for residue in chain:
                    if atom_selection in residue:
                        atoms.append(residue[atom_selection])
        return atoms

    def extract_chain_coords(self, struct, chain_id, atom_selection="CA"):
        atoms = []
        for model in struct:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if atom_selection in residue:
                            atoms.append(residue[atom_selection].get_coord())
        return np.array(atoms)

    def chain_rmsd(self, struct1, chain1, struct2, chain2):
        # Extract coordinates for chain A (or another chain) and atom selection (e.g., C-alpha)
        coords1 = self.extract_chain_coords(struct1, chain1)
        coords2 = self.extract_chain_coords(struct2, chain2)

        diff = coords1 - coords2
        return np.sqrt(np.sum(diff ** 2) / len(coords1))


    def align(self, temp_sym, full_structure):
        # move structure to the center
        gom = temp_sym.get_vrt("VRT0")._vrt_orig
        full_structure.transform(rot=np.identity(3), tran = -gom)

        # orient to global coordinates
        self.align_along_z(temp_sym, full_structure)
        self.align_along_x(full_structure)

        # if align has been parsed then we orient it so the flip along the x-axis matches the structure to align to.
        if self.align_struct:
            full_structure2 = full_structure.copy()
            self.flip_around_x(full_structure2)
            self.align_along_x(full_structure2)
            full_structure1_rmsd = self.chain_rmsd(full_structure, self.a, self.align_struct, self.align_chain)
            full_structure2_rmsd = self.chain_rmsd(full_structure2, self.a, self.align_struct, self.align_chain)
            if full_structure2_rmsd < full_structure1_rmsd:
                full_structure = full_structure2
        return full_structure

    def flip_around_x(self, full_structure):
        # flip it around the x axis
        rot = rotation_matrix([1, 0, 0], 180)
        full_structure.transform(rot=rot, tran=[0, 0, 0])

    def align_along_z(self, temp_sym, full_structure):
        # move the structure along the z direction
        current_z_axis = temp_sym.get_vrt("VRT0")._vrt_z
        # this (.transform()) right multiplies and CoordinateFrame.rotate() left multiplies so we have to be carefull.
        # This is the right order, for right multiplication as current_z_axis should go onto [0, 0, 1]
        rot = rvv( [0, 0, 1], current_z_axis)
        full_structure.transform(rot=rot, tran = [0, 0, 0])

    def align_along_x(self, full_structure):
        # move the center of mass of chain A so it points towards the x axis
        chain_a = full_structure[0][self.a]
        com = self.get_center_of_mass_of_chain(chain_a)
        rot = rvv( [1, 0, 0 ], com) # align
        full_structure.transform(rot=rot, tran = [0, 0, 0])

    def get_structure(self):
        # Use make_symmdef_file.pl to create a regular symmetry file
        temp_struct, temp_sym = self.make_regular_symmetry()
        parser = PDBParser()
        full_structure = parser.get_structure('Protein', temp_struct)

        full_structure = self.align(temp_sym, full_structure)

        io = PDBIO()
        io.set_structure(full_structure)
        # save full structure
        io.save(f"{self.out}/{Path(self.struct_file).stem}_norm_symm.pdb")
        # extract main pose
        with tempfile.NamedTemporaryFile(mode='w+') as f:
            io.save(f, ChainSelect(self.a))
            pose_a = pose_from_file(f.name)
        i_pose, b_pose = None, None
        # extract i pose (if defined)
        if self.i is not None:
            io.set_structure(full_structure)
            with tempfile.NamedTemporaryFile(mode='w+') as f:
                io.save(f, ChainSelect(self.i))
                i_pose = pose_from_file(f.name)
        # extract b pose (if defined)
        if self.b is not None:
            io.set_structure(full_structure)
            with tempfile.NamedTemporaryFile(mode='w+') as f:
                io.save(f, ChainSelect(self.b))
                b_pose = pose_from_file(f.name)
        return full_structure, pose_a, i_pose, b_pose

    def get_geometric_center(self, struct):
        xyzs = []
        for model in struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.get_fullname().replace(" ", "") == "CA":
                            xyzs.append(list(atom.coord))
        return sum(np.array(xyzs)) / len(xyzs)

    def extract_struct_coordinates(self, struct):
        # Extract coordinates into the chain_map
        chain_map = {}
        for model in struct:
            for chain in model:
                chain_map[chain.id] = []
                for residue in chain:
                    for atom in residue:
                        if atom.get_fullname().replace(" ", "") == "CA":
                            chain_map[chain.id].append(list(atom.coord))
        chain_map = {c: np.array(xyz) for c, xyz in chain_map.items()}
        return chain_map

    def get_anchor_residue(self, pose):
        return residue_center_of_mass(pose.conformation(), 1, pose.chain_end(1))

    def get_system_values(self, num_chains):
        if self.point_group == "C":
            if self.chains_type == "minimal":
                n_chains_to_include, multiplier = 2, num_chains
            if self.chains_type == "unique":
                n_chains_to_include = math.ceil(num_chains / 2)
                multiplier = num_chains
            if self.chains_type == "full":
                n_chains_to_include, multiplier = num_chains, 1
        else:
            raise NotImplementedError
        return n_chains_to_include, multiplier

    def make_cyclical_symmetry(self):
        symmetrical_struct, a_pose, i_pose, _ = self.get_structure()
        num_chains = len(list(symmetrical_struct.get_chains()))
        anchor_resi = self.get_anchor_residue(a_pose)

        # create the symmetrical setup
        setup = SymmetrySetup()

        # global center VRT
        center_vrt = CoordinateFrame("VRT_global", [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0])
        setup.add_vrt(center_vrt)

        # make a vrt that:
        # 1. Is level to the anchor CA atom (z coordinate is 0)
        # 2. Points with its x_vrt towards the the anchor atom (-x bug in Rosetta we rotate yet another 180 degrees)
        # 1:
        VRT1_x_t_ref = setup.copy_vrt("VRT_global", "VRT1_x_t_ref")
        anchor_CA_atom_pose_a = np.array(a_pose.residue(anchor_resi).atom("CA").xyz())
        trans = vector_projection(anchor_CA_atom_pose_a, center_vrt._vrt_z)
        VRT1_x_t_ref.translate(trans)
        assert np.isclose(VRT1_x_t_ref._vrt_z, center_vrt._vrt_z).all()
        # 2:
        onto_x_vec = normalize(anchor_CA_atom_pose_a - VRT1_x_t_ref._vrt_orig)
        onto_x_vec = - onto_x_vec
        rot = rvv(VRT1_x_t_ref._vrt_x , onto_x_vec)
        VRT1_x_t_ref.rotate(rot)
        assert np.isclose(VRT1_x_t_ref._vrt_z, center_vrt._vrt_z).all()
        setup.add_vrt(VRT1_x_t_ref)
        setup.add_jump("JUMP_global_to_1", "VRT_global", "VRT1_x_t_ref")


        # connect to the reference
        setup.add_vrt(setup.copy_vrt("VRT1_x_t_ref", "VRT1_x_t"))
        setup.add_jump("JUMP1_x_t", "VRT1_x_t_ref", "VRT1_x_t")

        # make 3 reference dofs and their connections
        setup.add_vrt(setup.copy_vrt("VRT1_x_t_ref", "VRT1_x_r_ref"))
        setup.add_vrt(setup.copy_vrt("VRT1_x_t_ref", "VRT1_x_r"))
        setup.add_jump("JUMP1_x_r_ref", "VRT1_x_t", "VRT1_x_r_ref")
        setup.add_jump("JUMP1_x_r", "VRT1_x_r_ref", "VRT1_x_r")

        setup.add_vrt(setup.copy_vrt("VRT1_x_t_ref", "VRT1_y_r_ref"))
        setup.add_vrt(setup.copy_vrt("VRT1_x_t_ref", "VRT1_y_r"))
        setup.add_jump("JUMP1_y_r_ref", "VRT1_x_r", "VRT1_y_r_ref")
        setup.add_jump("JUMP1_y_r", "VRT1_y_r_ref", "VRT1_y_r")

        setup.add_vrt(setup.copy_vrt("VRT1_x_t_ref", "VRT1_z_r_ref"))
        setup.add_vrt(setup.copy_vrt("VRT1_x_t_ref", "VRT1_z_r"))
        setup.add_jump("JUMP1_z_r_ref", "VRT1_y_r", "VRT1_z_r_ref")
        setup.add_jump("JUMP1_z_r", "VRT1_z_r_ref", "VRT1_z_r")

        # make an sds vrt just before the last jump
        setup.add_vrt(setup.copy_vrt("VRT1_z_r", "VRT1_sds"))
        setup.add_jump("JUMP1_to_sds", "VRT1_z_r", "VRT1_sds")

        # make final jump to the SUBUNIT
        setup.add_jump("JUMP1_to_subunit", "VRT1_sds", "SUBUNIT")

        # Now we copy all of that and rotate it around the global z axis to create the other SUBUNIT:
        anchor_CA_atom_pose_i = np.array(i_pose.residue(anchor_resi).atom("CA").xyz())
        angle = vector_angle(anchor_CA_atom_pose_a - VRT1_x_t_ref._vrt_orig, anchor_CA_atom_pose_i - VRT1_x_t_ref._vrt_orig)
        n_chains_to_include, multiplier = self.get_system_values(num_chains)
        vrt_copy = list(setup.vrts)
        jumps_copy = dict(setup.jumps)
        for chain in range(2, n_chains_to_include + 1):
            rot = rotation_matrix(center_vrt._vrt_z, angle * (chain -1) )
            for vrt in vrt_copy: # doing to deepcopy
                if vrt.name == "VRT_global":
                    continue
                vrt = setup.copy_vrt(vrt.name, vrt.name.replace("1", f"{chain}"))
                vrt.rotate(rot)
                setup.add_vrt(vrt)
            for jumpname, vrts in jumps_copy.items(): # doing deepcopy
                new_jumpname = jumpname.replace("1", f"{chain}")
                new_src = vrts[0].replace("1", f"{chain}")
                new_dest = vrts[1].replace("1", f"{chain}")
                setup.add_jump(new_jumpname, new_src, new_dest)

        # setup the 4 degrees of freedom
        setup.add_dof(f"JUMP1_x_t", 'x', "translation", np.linalg.norm(anchor_CA_atom_pose_a - VRT1_x_t_ref._vrt_orig))
        setup.add_dof(f"JUMP1_x_r", 'x', "rotation", 0)
        setup.add_dof(f"JUMP1_y_r", 'y', "rotation", 0)
        setup.add_dof(f"JUMP1_z_r", 'z', "rotation", 0)

        # n_chains_to_include = 1
        # setup the jumpgroups
        setup.add_jumpgroup("JUMPGROUP1", *(f"JUMP{i}_x_t" for i in range(1, n_chains_to_include + 1)))
        setup.add_jumpgroup("JUMPGROUP2", *(f"JUMP{i}_x_r" for i in range(1, n_chains_to_include + 1)))
        setup.add_jumpgroup("JUMPGROUP3", *(f"JUMP{i}_y_r" for i in range(1, n_chains_to_include + 1)))
        setup.add_jumpgroup("JUMPGROUP4", *(f"JUMP{i}_z_r" for i in range(1, n_chains_to_include + 1)))
        setup.add_jumpgroup("JUMPGROUP5", *(f"JUMP{i}_to_subunit" for i in range(1, n_chains_to_include + 1)))

        # set name, energy lines, anchor residue, and headers
        setup.symmetry_name = f"{Path(self.struct_file).stem}_ref_symm"
        setup.energies = f"{num_chains}*VRT1_sds + " + " + ".join([f"{multiplier}*(VRT1_sds:VRT{i}_sds)" for i in range(1, n_chains_to_include + 1)])
        setup.anchor = anchor_resi
        setup.headers["symmetry_type"] = f"C{num_chains}"
        setup.headers["actual_chains"] = num_chains
        setup.headers["chains_represented"] = n_chains_to_include

        # output the symmetry
        with open(f"{self.out}/{Path(self.struct_file).stem}_ref.symm", "w") as f:
            f.write(setup.make_symmetry_definition())

        # move the pose to origo and dump it
        setup.recenter = True
        a_pose.translate(-anchor_CA_atom_pose_a)
        a_pose.dump_pdb(f"{self.out}/{Path(self.struct_file).stem}_ref_INPUT.pdb")


    def make_regular_symmetry(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(f"{d}/input.pdb")
            shutil.copy(Path(self.struct_file).resolve(), p)
            if self.point_group == "C":
                command = f"cd {d} && {self.make_symmdef_file} -p {p} -a {self.a} -i {self.i} > symdef.symm"
            else:
                raise NotImplementedError
            subprocess.run(command, shell=True)
            # move struct    file
            symmetrical_struct = Path(f"{d}/{p.stem}_symm.pdb")
            new_symmetrical_struct = f"{self.out}/{symmetrical_struct.name}"
            shutil.move(symmetrical_struct, new_symmetrical_struct)
            # move symmetrical file
            symdef = f"{d}/symdef.symm"
            ss = SymmetrySetup(symdef=symdef)
        return new_symmetrical_struct, ss
