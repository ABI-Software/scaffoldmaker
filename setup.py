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
    "opencmiss.utils @ https://api.github.com/repos/OpenCMISS-Bindings/opencmiss.utils/tarball/master",
    "ZincPythonTools @ https://api.github.com/repos/OpenCMISS-Bindings/ZincPythonTools/tarball/master",
    "scipy",
    "numpy",
]
source_license = readfile("LICENSE")

setup(
    name="scaffoldmaker",
    version="0.1.2",
    description="Python client for generating anatomical scaffolds using OpenCMISS-Zinc",
    long_description="\n".join(readme) + source_license,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Medical Science Apps."
    ],
    author="Auckland Bioengineering Institute",
    author_email="r.christie@auckland.ac.nz",
    url="https://github.com/ABI-Software/scaffoldmaker",
    license="Apache Software License",
    packages=find_packages("src"),
    package_dir={"": "src"},
    include_package_data=True,
    zip_safe=False,
    install_requires=requires,
    )
