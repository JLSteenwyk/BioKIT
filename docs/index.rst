.. image:: _static/img/logo.png
   :width: 55%
   :align: center
   :target: https://jlsteenwyk.com/BioKIT

^^^^^


BioKIT, a versatile toolkit for processing or conducting analyses
on alignments, coding sequences, fastq files, and genome files.

If you found BioKIT useful, please cite *BioKIT: a versatile toolkit
for processing and analyzing diverse types of sequence data*. bioRxiv. doi: |doiLink|_.

.. _doiLink: https://www.biorxiv.org/content/10.1101/2021.10.02.462868v2
.. |doiLink| replace:: 10.1101/2021.10.02.462868

Quick Start
-----------
**1) Installation**

To install using *pip*, we strongly recommend building a virtual environment to avoid 
software dependency issues. To do so, execute the following commands:

.. code-block:: shell

	# create virtual environment
	python -m venv .venv
	# activate virtual environment
	source .venv/bin/activate
	# install biokit
	pip install jlsteenwyk-biokit

**Note, the virtual environment must be activated to use biokit.**

After using BioKIT, you may wish to deactivate your virtual environment and can do so using the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

|

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the
following commands:

.. code-block:: shell

	# download
	git clone https://github.com/JLSteenwyk/BioKIT.git
	cd BioKIT/
	# create virtual environment
	python -m venv .venv
	# activate virtual environment
	source .venv/bin/activate
	# install
	make install

To deactivate your virtual environment, use the following command:

.. code-block:: shell

	# deactivate virtual environment
	deactivate

**Note, the virtual environment must be activated to use biokit.**

|

To install via anaconda, execute the following command:

.. code-block:: shell

	conda install -c jlsteenwyk jlsteenwyk-biokit

Visit here for more information:
https://anaconda.org/JLSteenwyk/jlsteenwyk-biokit

|

**2) Usage**

Get the help message from BioKIT:

.. code-block:: shell

	biokit -h

|

^^^^


.. toctree::
	:maxdepth: 4

	about/index
	usage/index
	tutorials/index
	change_log/index
	other_software/index
	frequently_asked_questions/index

^^^^
