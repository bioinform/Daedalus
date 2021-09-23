
from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='extract-umi',
    version='0.1.0',
    packages=["extractUMI"],
    entry_points={
        'console_scripts': [
            'extract-umi = extractUMI.extractUMI:main'
        ],
    },
    url='http://ghe-rss.roche.com/pls-red-packages/extract-umi',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for extracting UMI information from read headers, adding UMI columns to an existing tsv file',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'pandas',
        'argparse'
        ]
)


