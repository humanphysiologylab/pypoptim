from setuptools import setup, find_packages

setup(
    name="pypoptim",
    version="0.12",
    packages=find_packages(),
    url="https://github.com/humanphysiologylab/pypoptim",
    author="Andrey Pikunov, Roman Syunyaev",
    author_email="pikunov@phystech.edu, roman.syunyaev@gmail.com",
    install_requires=[
        "numpy>=2.2.3,<3.0", 
        "pandas>=2.2.3,<3.0",
        "scikit-learn>=1.7.2,<1.8.0",
        "pytest>=8.3.5,<9.0",
        "numba>=0.62.1,<0.63.0"],
    python_requires='>=3.13,<3.14', 
)
