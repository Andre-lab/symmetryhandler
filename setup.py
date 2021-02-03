from distutils.core import setup

setup(
    name='symmetryhandler',
    version='0.1',
	scripts=['scripts/visualize_symmetry.py'],
    packages=['symmetryhandler'],
    # package_dir={'':'symmetryhandler'}
    url='https://github.com/Andre-lab/symmetryhandler',
    license='MIT',
    author='mads',
    author_email='mads.jeppesen@biochemistry.lu.se',
    description='Simple utility library to handcraft symmetry files in Rosetta and visualizing them in PyMOL',
	install_requires=['numpy']
)
