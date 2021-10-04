<p align="center">
  <a href="https://github.com/jlsteenwyk/biokit">
    <img src="https://github.com/JLSteenwyk/BioKIT/blob/main/docs/_static/img/biokit_logo.jpg" alt="Logo" width="400">
  </a>
  <p align="center">
    <a href="https://jlsteenwyk.com/BioKIT/">Docs</a>
    ·
    <a href="https://github.com/jlsteenwyk/biokit/issues">Report Bug</a>
    ·
    <a href="https://github.com/jlsteenwyk/biokit/issues">Request Feature</a>
  </p>
    <p align="center">
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <a href="https://pypi.org/project/jlsteenwyk-biokit/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/jlsteenwyk-biokit">
        </a>
        <a href="https://github.com/JLSteenwyk/BioKIT/actions" alt="Build">
            <img src="https://img.shields.io/github/workflow/status/JLSteenwyk/BioKIT/CI">
        </a>
        <a href="https://codecov.io/gh/JLSteenwyk/BioKIT" alt="Coverage">
          <img src="https://codecov.io/gh/JLSteenwyk/BioKIT/branch/main/graph/badge.svg?token=5X9C6YAVIG">
        </a>
        <a href="https://github.com/jlsteenwyk/biokit/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/biokit">
        </a>
        <a href="https://twitter.com/intent/follow?screen_name=jlsteenwyk" alt="Author Twitter">
            <img src="https://img.shields.io/twitter/follow/jlsteenwyk?style=social&logo=twitter"
                alt="follow on Twitter">
        </a>
    </p>
</p>

BioKIT is a UNIX shell toolkit for processing molecular sequence data.<br /><br />
If you found biokit useful, please cite the following: *BioKIT: a versatile toolkit for processing and analyzing diverse types of sequence data*. bioRxiv. doi: [10.1101/2021.10.02.462868](https://www.biorxiv.org/content/10.1101/2021.10.02.462868v1).
<br /><br />

---

This documentation covers downloading and installing BioKIT. Details about each function as well as tutorials for using BioKIT are available in the <a href="https://jlsteenwyk.com/BioKIT/">online documentation</a>.

<br />

**Installation** <br />

**If you are having trouble installing BioKIT, please contact the lead developer, Jacob L. Steenwyk, via [email](https://jlsteenwyk.com/contact.html) or [twitter](https://twitter.com/jlsteenwyk) to get help.**

To install using *pip*, we strongly recommend building a virtual environment to avoid software dependency issues. To do so, execute the following commands:
```shell
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install biokit
pip install jlsteenwyk-biokit
```

**Note, the virtual environment must be activated to use biokit.**

After using biokit, you may wish to deactivate your virtual environment and can do so using the following command:
```shell
# deactivate virtual environment
deactivate
```

<br />

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the following commands:
```shell
# download
git clone https://github.com/JLSteenwyk/BioKIT.git
cd biokit/
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install
make install
```
To deactivate your virtual environment, use the following command:
```shell
# deactivate virtual environment
deactivate
```
**Note, the virtual environment must be activated to use biokit.**

<br />

To install via anaconda, execute the following command:
```shell
conda install -c jlsteenwyk jlsteenwyk-biokit
```
Visit here for more information:
https://anaconda.org/JLSteenwyk/jlsteenwyk-biokit

<br />

To test biokit installation, launch the help message

```shell
biokit -h
```
