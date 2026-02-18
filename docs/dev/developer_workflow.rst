.. _developer_workflow:

Developer Workflow
==================

This page summarizes common local quality and profiling workflows.

Environment
-----------

.. code-block:: shell

   python -m venv venv
   source venv/bin/activate
   make dev.setup

Core Quality Checks
-------------------

.. code-block:: shell

   make lint
   make typecheck
   make guard.private-api
   make test.fast
   make quality

Integration Suite
-----------------

.. code-block:: shell

   make test.integration
   venv/bin/pytest -q tests/integration/test_output_format_matrix.py
   venv/bin/pytest -q tests/integration/test_cli_output_contracts.py

Benchmarking Hotspots
---------------------

Run benchmark summary in terminal:

.. code-block:: shell

   make bench.profile

Tune iterations/warmups:

.. code-block:: shell

   N=10 WARMUPS=2 make bench.profile

Save benchmark output for comparison:

.. code-block:: shell

   make bench.profile.save

Saved output path:

.. code-block:: text

   output/bench_profile.tsv

Notes
-----

- Benchmark values can vary by machine load; compare medians across repeated runs.
- The benchmark command currently tracks:
  - ``genome_assembly_metrics``
  - ``character_frequency``
  - ``relative_synonymous_codon_usage``
  - ``translate_sequence``
