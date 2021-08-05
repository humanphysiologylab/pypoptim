from setuptools import setup

setup(
    name="pypoptim",
    version="0.11",
    packages=["pypoptim"],
    url="https://github.com/humanphysiologylab/pypoptim",
    author="Andrey Pikunov",
    author_email="pikunov@phystech.edu",
    install_requires=["numpy", "pandas", "scikit-learn", "numba", "pytest"],
)
