Overview
=================================================
geostochpy is a project that was built with the need for
Being able to generate stochastic scenarios of earthquakes based on geometry in subduction areas
along the trench.There has been an approach in the Chilean subduction zone.The project includes functions
To generate grids, to generate SLIP distributions with different TAPERS, to obtain deformations, and for the
Pre-processing of data for tsunami softwares as a comcot

How does work?
-------------------------------------------------------
Actually, the Stochastic Slip generation its based on Karhunen-Loéve expansion

.. math::

    X_t=\sum_{k=1}^\infty Z_ke_k(t) 

For the stochastic generation of slips:

.. math:: 

    Slip=\mu+\sum_{k=1}^Nz_k\sqrt{\lambda_k}v_k

That a stochastic process can be represented as an infinite linear combination of orthogonal functions.

In contrast to a Fourier series :math:`F(t) = a_0 + \sum_{n=1}^{\infty} \left( a_n \cos\left(\frac{2\pi n t}{T}\right) + b_n \sin\left(\frac{2\pi n t}{T}\right) \right)`
where the coefficients are fixed numbers and the expansion basis consists
of sinusoidal functions (that is, sine and cosine functions), the coefficients in the Karhunen–Loève 
theorem are random variables and the expansion basis depends on the process.
In fact, the orthogonal basis functions used in this representation are determined by the covariance function of the process.

Well, the stochastic generation is based on the Paper of :cite:t:`LeVeque2017`
Where the covariance matrix is

.. math:: 

    C_{ij}=exp(-(d_{strike}(i,j)/r_{strike})-(d_{dip}(i,j)/r_{dip}))

Where

.. math:: d_{Strike} (i, j), d_{dip}(i, j)

are estimated of the distance between subfaults i and j in the strike and dip respectively, and  

.. math:: r_{Strike}, r_{Dip}

are the lengths of correlation in each direction.

It is defined :math:`d_{dip}(i, j)`
using the in depth difference between two subfaults and the angle of Dip as
:math:`d_{dip}(i,j)=d_{depth}/sin(dip)`
.Setting :math:`d_{strike}=\sqrt{d_{ij}^2-d_{dip}(i,j)^2}`

Also, for the generation was implemented an 2-D filter based on Tukey Window along the strike and along the dip.
For more details you can see the :ref:`gallery-section`

.. bibliography::
