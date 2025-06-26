Working on documentation
=====================================================================

Instruction for local doc build after you changes rst files 
-----------------------------------------------------------

0. If you have not created sphinx python environment, please create it first 

::

    conda env create -f sphinx.yaml 

1. Activate sphinx environment

::

    conda activate sphinx

2. in your docs directory run:

::

    make html
    python3 -m http.server

3. Open your browser and copy ``http://localhost:8000/`` in address bar to see the docs in html

