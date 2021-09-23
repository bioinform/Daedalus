from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='bam-streamer',
    version='0.1.0',
    packages=["BamStreamer"],
    url='http://ghe-rss.roche.com/pls-red-packages/bam-streamer',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A package for reading one or more bam Files iteraing fastq file records for both single-end or paired-end fastq files.',
    long_description=long_description,
    long_description_content_type="text/markdown"
)
