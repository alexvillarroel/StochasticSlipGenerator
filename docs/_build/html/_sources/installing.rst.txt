Installing
==================================================================================================
For install this tool, it's so simple, you need mamba or conda package manager, which takes care of setting up a virtual environment.

First you need to clone the repository

.. _repo:
    
    git clone https://github.com/alexvillarroel/StochasticSlipGenerator

Once downloaded, get into folder and write:

.. tab-set::

    .. tab-item:: mamba
        :sync: mamba

        ::

            mamba env create -f environment.yml

    .. tab-item:: conda
        :sync: conda

        ::

            conda env create -f environment.yml

Now you can make a test using an example of :ref:`gallery-section`