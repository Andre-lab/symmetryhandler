from scipy.optimize import dual_annealing
import numpy as np
from Bio.PDB import Superimposer
from symmetryhandler.mathfunctions import rotation_matrix
import copy
from scipy.optimize import minimize

class SymmetryOptimizer:

    def __init__(self, target, reference, target_chain="A", reference_chain="A", translation_bounds=(-100, 100), rotation_bounds=(-180, 180)):
        self.target = copy.deepcopy(target)
        self.reference = reference
        self.target_chain = target_chain
        self.reference_chain = reference_chain
        self.bounds = (translation_bounds, rotation_bounds)

    def extract_chain_coords(self, struct, chain_id, atom_selection="CA"):
        atoms = []
        for model in struct:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if atom_selection in residue:
                            atoms.append(residue[atom_selection].get_coord())
        return np.array(atoms)

    def optimize(self):
        # Output the start RMSD
        # print("Starting RMSD:", self.rmsd(self.target_coords, self.reference_coords))

        # Run the optimization
        best_params = dual_annealing(self.compute_rmsd, self.bounds)

        # Output the Final rmsd
        best_rmsd =self.compute_rmsd(best_params.x)
        # print("Final RMSD:", best_rmsd)

        return best_params, best_rmsd

    def apply(self):
        """Runs dual annealing to minimize RMSD."""
        # Create starting coordinates
        self.target_coords = self.extract_chain_coords(self.target, self.target_chain)
        self.reference_coords = self.extract_chain_coords(self.reference, self.reference_chain)

        # Run the same optimization
        best_params, rmsd = self.optimize()

        # Flip the coordinates
        self.target_coords = np.dot(self.target_coords, rotation_matrix([1, 0, 0], 180))

        # Run the same optimization
        best_params_flip, rmsd_flip = self.optimize()

        if rmsd <= rmsd_flip:
            return self.apply_best_params(best_params), rmsd
        else:
            self.flip_structure(self.target)
            return self.apply_best_params(best_params_flip), rmsd_flip

    def flip_structure(self, structure):
        rot = rotation_matrix([1, 0, 0], 180)
        structure.transform(rot=rot, tran=[0,0,0])

        # potential final minimize
        # result_lbfgsb = minimize(
        #     self.compute_rmsd,  # Objective function
        #     best_params.x,  # Initial guess from dual_annealing
        #     method="L-BFGS-B",
        #     bounds=bounds
        # )
        #
        # rmsd = self.compute_rmsd(result_lbfgsb.x)
        # print("Final RMSD:", rmsd)



    def apply_best_params(self, best_params):
        transformed_structure = copy.deepcopy(self.target)
        translation = best_params.x[0]
        angle = best_params.x[1]
        rot = rotation_matrix([0, 0, 1], angle)
        transformed_structure.transform(rot=rot, tran=[0, 0, translation])
        return transformed_structure

    def x_flip(self, structure, params):
        # turn the continuous 0-1 variable into True/False (1/0)
        flip = int(round(params[2]))
        if flip:
            rot = rotation_matrix([1, 0, 0], 180)
            structure.transform(rot=rot, tran=[0, 0, 0])
        return structure

    @staticmethod
    def z_translation(coords, translation):
        """Translate structure by shift vector."""
        coords[:, 2] += translation
        return coords

    @staticmethod
    def z_rotation(coords, angle):
        """Rotate structure around a random axis by a given angle."""
        rot = rotation_matrix([0, 0, 1], angle)
        return np.dot(coords, rot)

    def rmsd(self, coords1, coords2):
        # Extract coordinates for chain A (or another chain) and atom selection (e.g., C-alpha)
        diff = coords1 - coords2
        return np.sqrt(np.sum(diff ** 2) / len(coords1))

    # todo optimization: Keep everything as a pointcloud and just move them. Should be fairly easy to do.
    def compute_rmsd(self, params):
        """Objective function: Compute RMSD after applying a random movement."""
        coords = copy.deepcopy(self.target_coords)
        coords = self.z_translation(coords, params[0])
        coords = self.z_rotation(coords, params[1])
        rmsd = self.rmsd(coords, self.reference_coords)
        return rmsd
