from os import path
from setuptools import setup, find_packages

from biokit.version import __version__

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as f:
    long_description = f.read()

CLASSIFIERS = [
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering",
]

REQUIRES = ["biopython==1.76", "numpy==1.18.2", "cython"]

setup(
    name="biokit",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/biokit",
    packages=find_packages(),
    classifiers=CLASSIFIERS,
    entry_points={
        "console_scripts": [
            "biokit = biokit.biokit:main",
            # "bk_aln_len = biokit.biokit:alignment_length", # Alignment-based functions
        ]
    },
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
)

# push new version to pypi
# rm -rf dist
# python setup.py sdist bdist_wheel --universal
# twine upload dist/* -r pypi
# then push to anaconda
