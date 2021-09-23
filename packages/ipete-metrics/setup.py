from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='ipete-metrics',
    version='0.1.2',
    packages=["ipeteMetrics"],
    entry_points={
        'console_scripts': [
            'ipete-metrics = ipeteMetrics.ipeteMetrics:main'
        ],
    },
    url='http://ghe-rss.roche.com/pls-red-packages/ipete-metrics',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for defining Diversity and Entropy Metrics from ImmunoPETE libraries ',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'pandas',
        'argparse',
        'numpy'
        ]
)
