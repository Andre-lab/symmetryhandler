from shapedesign.src.utilities.tests import setup_test
from symmetryhandler.symmetryhandler import SymmetrySetup
from shapedesign.src.utilities.kinematics import perturb_jumpdof
from shapedesign.src.utilities.pose import get_position_info
import time
# from shapedesign.src.utilities.pose import set_all_dofs_to_zero
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit
from pyrosetta import Pose



def test_visualization():
    setup = SymmetrySetup()
    setup.read_from_file("inputs/1stm.symm")
    setup.visualize()
    setup.print_visualization("out/1stm.symm")

def test_read_from_pose():

    pose, pmm, cmd, symm_file = setup_test(return_symmetry_file=True)
    pmm.keep_history(True)
    pmm.apply(pose)

    position_info = get_position_info(pose, dictionary=True)

    setup = SymmetrySetup()
    setup.read_from_file(symm_file)

    perturb_jumpdof(pose, "JUMP5fold1", 3, 10)
    perturb_jumpdof(pose, "JUMP5fold111", 1, 10)
    # set_jumpdof(pose, "JUMP5fold1111_subunit", 6, -10)

    setup.update_from_pose(pose)

    pmm.apply(pose)
    setup.visualize(ip="10.8.0.14", port="9123")

def set_all_dofs_to_zero(pose):
    from shapedesign.src.utilities.kinematics import set_jumpdof
    from shapedesign.src.utilities.pose import  dof_map
    for jump, dofs in get_position_info(pose, dictionary=True).items():
        # if jump == "JUMP5fold1111_subunit":
        #     continue
        for dof, old_val in dofs.items():
            set_jumpdof(pose, jump, dof_map[dof], 0)

def test_create_independent_icosahedral_symmetries():
    pose_non_sym = setup_test("5cvz", symmetry=False, pymol=False)
    pose, pmm, cmd, symm_file = setup_test("5cvz", return_symmetry_file=True)
    pmm.keep_history(True)
    setup = SymmetrySetup()
    setup.read_from_file(symm_file)

    perturb_jumpdof(pose, "JUMP5fold1", 3, 10)
    perturb_jumpdof(pose, "JUMP5fold1", 6, 20)
    perturb_jumpdof(pose, "JUMP5fold1111", 4, 10)
    perturb_jumpdof(pose, "JUMP5fold1111", 5, 20)
    perturb_jumpdof(pose, "JUMP5fold1111", 6, 30)
    # # set_jumpdof(pose, "JUMP5fold1111_subunit", 6, -10)
    pmm.apply(pose)

    fold5_setup, fold3_setup, fold2_1_setup, fold2_2_setup = setup.create_independent_icosahedral_symmetries(pose)

    pose_dofs_are_0 = pose.clone()
    set_all_dofs_to_zero(pose_dofs_are_0)
    asymmetric_pose = Pose()
    extract_asymmetric_unit(pose_dofs_are_0, asymmetric_pose, False)

    #pmm.apply(pose_non_sym); pmm.apply(asymmetric_pose) #ARE IDENTICAL before symmetry = GOOD!

    # fold5
    # fold5_setup.visualize(ip="10.8.0.10", port="9123")
    fold5_setup.make_symmetric_pose(asymmetric_pose)
    pmm.apply(asymmetric_pose)

    # fold3
    # pose, pmm, cmd = setup_test("4v4m", symmetry=False)

    # fold3_setup.visualize(ip="10.8.0.10", port="9123")
    asymmetric_pose = Pose()
    extract_asymmetric_unit(pose_dofs_are_0, asymmetric_pose, False)
    fold3_setup.make_symmetric_pose(asymmetric_pose)
    pmm.apply(asymmetric_pose)

    # fold2_1
    # pose, pmm, cmd = setup_test("4v4m", symmetry=False)
    # fold2_1_setup.visualize(ip="10.8.0.10", port="9123")
    asymmetric_pose = Pose()
    extract_asymmetric_unit(pose_dofs_are_0, asymmetric_pose, False)
    fold2_1_setup.make_symmetric_pose(asymmetric_pose)
    pmm.apply(asymmetric_pose)

    # fold2_2
    # pose, pmm, cmd = setup_test("4v4m", symmetry=False)
    # fold2_2_setup.visualize(ip="10.8.0.10", port="9123")
    asymmetric_pose = Pose()
    extract_asymmetric_unit(pose_dofs_are_0, asymmetric_pose, False)
    fold2_2_setup.make_symmetric_pose(asymmetric_pose)
    pmm.apply(asymmetric_pose)

def test_rotations_to_2folds():
    pose, pmm, cmd, symm_file = setup_test("1stm", return_symmetry_file=True)
    setup = SymmetrySetup()
    setup.read_from_file(symm_file)
    setup.visualize(ip="10.8.0.10", port="9123")
    pmm.apply(pose)
    start = time.time()
    #a, b = setup.rotations_to_2folds(pose, visualize=True, cmd=cmd)
    a, b = setup.rotations_to_2folds(pose)
    print(time.time() - start)
    pass






    #setup.visualize(ip="10.8.0.14", port="9123")
