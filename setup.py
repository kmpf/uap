import os
from distutils.core import setup
from __future__ import print_function


# List of lists with tools to check for
check_tools = {
    "virtualenv": ["virtualenv", "--version"],
    "gcc": ["gcc", "--version"],
    "git": ["git", "--version"]
}

for tool in keys(check_tools):
    try:
        subprocess.call(check_tools[tool])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print("%s is required but it's not installed." % tool,
                  file=sys.stderr)
        else:
            # Something else went wrong while trying to run `wget`
            raise

# Create new Python virtualenv
subprocess.call(["virtualenv", "python_env"])
# Activate the created virtualenv
uap_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/python_env/bin/activate_this.py' % uap_path

requires = [
    'pyyaml', 
    'numpy', 
    'biopython', 
    'psutil'
]

setup(
    name='uap',
    version='0.0.1',
    author='Christoph Kaempf',
    author_email='christoph.kaempf@ufz.de',
    description='Universal Analysis Pipeline',
    # README.txt has to be in reStructuredText format
    long_description=open("README.txt").read(),
    platforms='',
    license='',
    url='https://github.com/kmpf/uap',
    scripts=['uap'],
    install_requires=requires
)
