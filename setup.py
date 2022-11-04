import io
import os

from setuptools import setup, find_packages

SETUP_DIR = os.path.dirname(os.path.abspath(__file__))


def readfile(filename, split=False):
    with io.open(filename, encoding="utf-8") as stream:
        if split:
            return stream.read().split("\n")
        return stream.read()


readme = readfile("README.rst", split=True)
readme.append('License')
readme.append('=======')
readme.append('')
readme.append('::')
readme.append('')
readme.append('')

# For requirements not hosted on PyPi place listings
# into the 'requirements.txt' file.
requires = [
    # minimal requirements listing
    "opencmiss.maths",
    "opencmiss.utils >= 0.3",
    "opencmiss.zinc >= 3.9",
    "scipy",
    "numpy",
]
source_license = readfile("LICENSE")

setup(
    name="scaffoldmaker",
    version="0.8.0",
    description="Python client for generating anatomical scaffolds using OpenCMISS-Zinc",
    long_description="\n".join(readme) + source_license,
    long_description_content_type="text/x-rst",
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
