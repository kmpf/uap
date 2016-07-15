import os
import subprocess
from distutils.core import setup

# List of lists with tools to check for
# J: why checking only for these 3 tools?
check_tools = {
    "virtualenv": ["virtualenv", "--version"],
    "gcc": ["gcc", "--version"],
    "git": ["git", "--version"]
}

for tool in check_tools.keys():
    try:
        subprocess.call(check_tools[tool])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            sys.stderr.write("%s is required but it's not installed." % tool)
            sys.stderr.flush()
        else:
            # Something else went wrong while trying to run `wget`
            raise

# Create new Python virtualenv named 'python_env' 
# J: this creates a subdirectory 'python_env' containing all the python executable 
# J: files and libraries in a [bin,include,lib] style
subprocess.call(["virtualenv", "python_env"])

# Activate the created virtualenv
# J: is it deactivated at the end? is this neccessary?
uap_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/python_env/bin/activate_this.py' % uap_path

execfile(activate_this_file, dict(__file__=activate_this_file))

# J: here we require additional python libraries?
# J: are they checked 

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
    long_description=open("README.rst").read(),
    platforms='',
    license='',
    url='https://github.com/kmpf/uap',
    scripts=['uap'],
    install_requires=requires
)
