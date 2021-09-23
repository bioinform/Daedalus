from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='ipete-spikein-parser',
    version='0.1.0',
    packages=["ipeteSpikeins"],
    entry_points={
        'console_scripts': [
            'ipeteSpikeins = ipeteSpikeins.ipeteSpikeins:main'
        ],
    },
    url='http://ghe-rss.roche.com/pls-red-packages/ipete-spikein-parser',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for parsing immunoPETE spikein alignments, either splitting spikeins from normal or counting spikein fragments',
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
