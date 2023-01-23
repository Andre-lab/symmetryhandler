
def test_vrts_overlap_with_pose():
    from simpletestlib.test import setup_test
    from cubicsym.kinematics import randomize_all_dofs
    from cubicsym.cubicsetup import CubicSetup
    #for pdb in ("1STM", "6SS4"):
    #for pdb in ("1AEW", "1P3Y"):
    for pdb in ("7Q03", "1H0S"):
        pose, pmm, cmd, symdef = setup_test(name="T", file=pdb, return_symmetry_file=True, mute=True)
        cs = CubicSetup(symdef)
        for i in range(10):
            randomize_all_dofs(pose)
            cs.vrts_overlap_with_pose(pose, update_and_apply_dofs=True)

            # pmm.apply(pose)
            # cs.visualize(ip="10.8.0.10")

def test_apply_dofs():
    """Tests that the dofs follow the pose"""
    from simpletestlib.test import setup_test
    from cubicsym.kinematics import randomize_all_dofs
    from cubicsym.cubicsetup import CubicSetup
    #for pdb in ("1STM", "6SS4"):
    #for pdb in ("1AEW", "1P3Y"):
    for pdb in ("7Q03", "1H0S"):
        pose, pmm, cmd, symdef = setup_test(name="T", file=pdb, return_symmetry_file=True, mute=True)
        cs = CubicSetup(symdef)
        for i in range(10):
            randomize_all_dofs(pose)
            cs.vrts_overlap_with_pose(pose, update_and_apply_dofs=True)
            # cs.update_dofs_from_pose(pose)
            # pmm.apply(pose)
            # cs.visualize(ip="10.8.0.6")
            ...

# def test_visualization():
#     setup = SymmetrySetup()
#     setup.read_from_file("inputs/1stm.symm")
#     setup.visualize()
#     setup.print_visualization("out/1stm.symm")
#
# def test_read_from_pose():
#     pose, pmm, cmd, symm_file = setup_test(return_symmetry_file=True)
#     pmm.keep_history(True)
#     pmm.apply(pose)
#
#     position_info = get_position_info(pose, dictionary=True)
#
#     setup = SymmetrySetup()
#     setup.read_from_file(symm_file)
#
#     perturb_jumpdof(pose, "JUMP5fold1", 3, 10)
#     perturb_jumpdof(pose, "JUMP5fold111", 1, 10)
#     # set_jumpdof(pose, "JUMP5fold1111_subunit", 6, -10)
#
#     setup.update_dofs_from_pose(pose)
#
#     pmm.apply(pose)
#     setup.visualize(ip="10.8.0.14", port="9123")
#
#
# def test_output_vrts_as_pdb():
#     pose, pmm, cmd, symm_file = setup_test("1stm", return_symmetry_file=True)
#     setup = SymmetrySetup()
#     setup.output_vrts_as_pdb(pose)
#
#
# def test_add_angle_vrt():
#     pose, pmm, cmd, symm_file = setup_test("1stm", return_symmetry_file=True)
#     setup = SymmetrySetup()
#     setup.read_from_file(symm_file)
#     setup.add_angle_vrt(pose)
#
# def test_reference_symmetry():
#     from simpletestlib.test import setup_test
#     from symmetryhandler.symmetrysetup import SymmetrySetup
#     pose, pmm, cmd, symdef = setup_test(name="I", file="1STM", return_symmetry_file=True, mute=False)
#     assert SymmetrySetup(symdef).is_reference_symmetry()
#
