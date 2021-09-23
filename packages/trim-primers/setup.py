from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='trim-primers',
    version='0.1.0',
    packages=["trimPrimers"],
    entry_points={
        'console_scripts': [
            'trim-primers = trimPrimers.trimPrimers:main'
        ],
    },
    url='http://ghe-rss.roche.com/pls-red-packages/trim-primers',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for triming V and J primer alignments from read fragments.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'pandas',
        'pysam',
        'argparse',
        'fastq-streamer',
        'bam-streamer'
        ],
    dependency_links = [
        'git+http://ghe-rss.roche.com/pls-red-packages/fastq-streamer#egg=fastq-streamer-0.1.0',
        'git+http://ghe-rss.roche.com/pls-red-packages/bam-streamer#egg=bam-streamer-0.1.0'
    ] 
)
