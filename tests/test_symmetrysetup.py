def test_convert_to_reference_symmetry():
    from symmetryhandler.symmetrysetup import SymmetrySetup
    from pyrosetta import pose_from_file, init
    from simpletestlib.setup import setup_pymol
    init("-symmetry:initialize_rigid_body_dofs true -pdb_comments")

    pmm = setup_pymol()

    # # Cyclical symmetry
    # symdef = "inputs/cyclical/c10_single.symm"
    # symdef = "inputs/cyclical/c10.symm"
    symdef = "/home/mads/projects/production_symevodock/symmetrizing_tests/cyclical/c10.symm"
    ss = SymmetrySetup(symdef=symdef)
    # ss.convert_to_reference_symmetry()
    ss.visualize(ip="10.8.0.14")
    #pose = pose_from_file("inputs/cyclical/input_INPUT.pdb")
    pose = pose_from_file("/home/mads/projects/production_symevodock/symmetrizing_tests/cyclical/input_INPUT.pdb")
    ss.make_symmetric_pose(pose)
    pmm.apply(pose)

    # Dihedral symmetry
    symdef = "inputs/dihedral/d5.symm"
    ss = SymmetrySetup(symdef=symdef)
    # ss.convert_to_reference_symmetry()
    ss.visualize(ip="10.8.0.14")
    pose = pose_from_file("inputs/dihedral/input_INPUT.pdb")
    ss.make_symmetric_pose(pose)
    pmm.apply(pose)

    print("adf")