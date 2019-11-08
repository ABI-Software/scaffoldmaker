from setuptools import setup, find_packages
from setuptools.command.install import install
import os
import io

SETUP_DIR = os.path.dirname(os.path.abspath(__file__))

# List all of your Python package dependencies in the
# requirements.txt file

def readfile(filename, split=False):
    with io.open(filename, encoding="utf-8") as stream:
        if split:
            return stream.read().split("\n")
        return stream.read()

readme = readfile("README.rst", split=True)[3:]  # skip title
# For requirements not hosted on PyPi place listings
# into the 'requirements.txt' file.
requires = [
    # minimal requirements listing
    'opencmiss.utils @ https://api.github.com/repos/OpenCMISS-Bindings/opencmiss.utils/tarball/master',
    'ZincPythonTools @ https://api.github.com/repos/OpenCMISS-Bindings/ZincPythonTools/tarball/master',
    'scipy',
    'numpy',
]
source_license = readfile("LICENSE")

setup(name='scaffoldmaker',
    version='0.1.2',
    description='',
    long_description='\n'.join(readme) + source_license,
    classifiers=[
      "Development Status :: 3 - Alpha",
      "License :: OSI Approved :: Apache Software License",
      "Programming Language :: Python",
    ],
    author='Richard Christie',
    author_email='',
    url='',
    license='APACHE',
    packages=find_packages(exclude=['ez_setup',]),
    package_dir={'': 'src'},
    include_package_data=True,
    zip_safe=False,
    install_requires=requires,
    )
