from setuptools import setup, find_packages

setup(
    name="speczfind",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "astropy",
        "matplotlib"
    ],
    entry_points={
        "console_scripts": [
            "speczfind=speczfind.cli:main",
        ],
    },
)
