from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='Daedalus',
    version='0.1.0',
    url='http://ghe-rss.roche.com/plsRED-Bioinformatics/Daedalus',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A pipeline for Immune Repertoire sequencing analysis.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'pandas',
        'numpy',
        'sample-sheet==0.1.0'
        ]
)
