from symmetryhandler import SymmetrySetup
import argparse

parser = argparse.ArgumentParser("xml-rpc server")
parser.add_argument("file")
parser.add_argument("--apply_dofs", help="Will apply the dof movements as specified in the symmetry_file "
                                            "(same as doing '-symmetry:initialize_rigid_body_dofs true') in Rosetta"
                                         ".", default=True)
parser.add_argument("--mark_jumps", help="Will mark the jumps.", default=True)
parser.add_argument("--create_script", help="Creates a script with that name can be used to visualize the symmetry in PyMOL.",
                    default="")
parser.add_argument("--ip", default="localhost", help="Sets the ip address that the PyMOL server is at. ")
parser.add_argument("--port", defualt="9123", help="Sets the port the PyMOL server listens to.")
args = parser.parse_args()

if __name__ == "__main__":
    setup = SymmetrySetup()
    setup.read_from_file(args.file)
    if args.create_script:
        setup.print_visualization(args.create_script, args.apply_dofs)
    else:
        setup.visualize(args.apply_dofs, args.mark_jumps, args.ip, args.port)
