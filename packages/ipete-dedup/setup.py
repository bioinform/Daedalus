from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='ipete-dedup',
    version='0.1.1',
    packages=["ipeteDedup"],
    entry_points={
        'console_scripts': [
            'ipete-dedup = ipeteDedup.ipeteDedup:main'
        ],
    },
    url='http://ghe-rss.roche.com/pls-red-packages/ipete-dedup',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for deduplicating CDR3 sequences and defining quality based consensus in umi families',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'pandas',
        'argparse',
        'scipy',
        'SeqNetwork',
        'pysam'
        ],
    dependency_links = [
        'git+http://ghe-rss.roche.com/pls-red-packages/SeqNetwork#egg=SeqNetwork-0.1.0'
    ] 
)
