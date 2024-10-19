import random
from symmetryhandler.reference_kinematics import get_dofs, perturb_jumpdof_str_str

def randomize_dofs(pose, size):
    """Randomizes all the dofs in the given range of +/- the passed value. It picks a value uniformly."""
    for jump, dofs in get_dofs(pose).items():
        for dof, current_val in dofs.items():
            val = random.uniform(-size, +size)
            perturb_jumpdof_str_str(pose, jump, dof, val)
