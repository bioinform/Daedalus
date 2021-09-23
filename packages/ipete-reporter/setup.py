from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setup(
    name='ipete-reporter',
    version='0.1.2',
    packages=["ipeteReporter"],
    entry_points={
        'console_scripts': [
            'ipete-reporter = ipeteReporter.ipeteReporter:main'
        ],
    },
    url='http://ghe-rss.roche.com/pls-red-packages/ipete-reporter',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for Reporter Clone and Diversity Stats from ImmunoPETE pipeline runs',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'pandas',
        'SeqNetwork',
        'ipete-metrics',
        ],
    dependency_links = [
        'git+http://ghe-rss.roche.com/pls-red-packages/SeqNetwork#egg=SeqNetwork-0.1.0'
        'git+http://ghe-rss.roche.com/pls-red-packages/ipete-metrics#egg=ipete-metrics-0.1.2'
    ] 

)
