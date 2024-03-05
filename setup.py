from setuptools import setup, find_packages

setup(
    name='geostochpy',
    version='0.1',
    description='Geostochpy is a project that was built with the need for Being able to generate stochastic scenarios of earthquakes based on geometry in subduction areas along the trench.There has been an approach in the Chilean subduction zone.The project includes functions To generate grids, to generate SLIP distributions with different TAPERS, to obtain deformations, and for the Pre-processing of data for tsunami softwares as a comcot',
    author='Alex Villarroel Carrasco',
    author_email='alexvillarroel.ca@gmail.com',
    url='https://github.com/alexvillarroel/StochasticSlipGenerator',
    packages=find_packages(),
    package_data={
        'geostochpy': ['data/*'],
    }
)
