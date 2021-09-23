from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='parse-umi',
    version='0.1.0',
    packages=["parseUMI"],
    entry_points={
        'console_scripts': [
            'parse-umi = parseUMI.parseUMI:main'
        ],
    },
    url='http://ghe-rss.roche.com/pls-red-packages/parse-umi',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for labeling read pair headers with UMI information',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'pandas',
        'pysam',
        'argparse',
        'fastq-streamer',
        ],
    dependency_links = [
        'git+http://ghe-rss.roche.com/pls-red-packages/fastq-streamer#egg=fastq-streamer-0.1.0'
    ] 
)
