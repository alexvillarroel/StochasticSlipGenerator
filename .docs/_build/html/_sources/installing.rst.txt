Installing
==================================================================================================
For install this tool, it's so simple, you need mamba or conda package manager (mamba its my favorite), which takes care of setting up a virtual environment.

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

Next step its activate the environment with the following command:

.. tab-set::

    .. tab-item:: mamba
        :sync: mamba

        ::

            mamba activate geostochpy

    .. tab-item:: conda
        :sync: conda

        ::

            conda activate geostochpy

Once the environment is activated, you need to install the package, for this you need to go to the folder where the setup.py file is located and write:

.. code-block:: bash

    python setup.py install

Now you can make a test using an example of :ref:`gallery-section`