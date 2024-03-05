<<<<<<< HEAD
Geostochpy: A tool to create,filter and make plots of stochastic slip grids applied in Chile 🔧📚🌏🇨🇱
=======
geostochpy: A tool to create,filter and make plots of stochastic slip grids applied in Chile 🔧📚🌏🇨🇱
>>>>>>> 457a92b (Final Pre-Thesis Work)
=================================================
Overview
=================================================
geostochpy is a project that was built with the need for
Being able to generate stochastic scenarios of earthquakes based on geometry in subduction areas
along the trench.There has been an approach in the Chilean subduction zone.The project includes functions
To generate grids, to generate SLIP distributions with different TAPERS, to obtain deformations, and for the
Pre-processing of data for tsunami softwares as a comcot

How does work?
-------------------------------------------------------
Actually, the Stochastic Slip generation its based on Karhunen-Loéve expansion,
that explain a stochastic process can be represented as an infinite linear combination of orthogonal functions.

Installing
==================================================================================================
For install this tool, it's so simple, you need mamba or conda package manager, which takes care of setting up a virtual environment.

First you need to clone the repository

.. _repo:
    
    git clone https://github.com/alexvillarroel/StochasticSlipGenerator

Once downloaded, get into folder and write:

.. _repo:

    mamba env create -f environment.yml

or

.. _repo:

    conda env create -f environment.yml

Once you had created the environment, you need to activate it.

<<<<<<< HEAD
.. _repo:

    conda activate stochpy

And write 

.. _repo:

    python setup.py install

And now you can use the package!

For more information, you can visit the documentation in Stochpy Documentation https://alexvillarroel.github.io/StochasticSlipGenerator
=======
For more information, you can visit the documentation in :ref:` geostochpy Documentation <>`
>>>>>>> 457a92b (Final Pre-Thesis Work)


