from distutils.core import setup

setup(
    name='symmetryhandler',
    version='0.1',
	scripts=['symmetryhandler/scripts/visualize_symmetry.py'],
    packages=['symmetryhandler', 'symmetryhandler/scripts'],
    # package_dir={'':'symmetryhandler'}
    url='https://github.com/Andre-lab/symmetryhandler',
    license='MIT',
    author='mads',
    author_email='mads.jeppesen@biochemistry.lu.se',
    description='Simple utility library to handcraft symmetry files in Rosetta and visualizing them in PyMOL',
    package_data={
        'symmetryhandler/scripts': ['make_symmdef_file.pl'],
    },
	install_requires=['numpy']
)
