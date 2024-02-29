from setuptools import setup, find_packages

setup(
    name="pypoptim",
    version="0.11",
    packages=find_packages(),
    url="https://github.com/humanphysiologylab/pypoptim",
    author="Andrey Pikunov",
    author_email="pikunov@phystech.edu",
    install_requires=["numpy", "pandas", "scikit-learn==1.3.2", "pytest", "numba"],
)
