"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyigrf",
    version="0.0.1",
    author="Pablo Reyes",
    author_email="pabloreyesfirpo@gmail.com",
    description="IGRF",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pmreyes2/pyigrf",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
