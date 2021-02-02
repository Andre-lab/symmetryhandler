from symmetryhandler import SymmetrySetup

def test_visualization():
    setup = SymmetrySetup()
    setup.read_from_file("inputs/1stm.symm")
    setup.visualize()
    setup.print_visualization("out/1stm.symm")

