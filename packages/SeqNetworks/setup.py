from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='SeqNetwork',
    version='0.1.0',
    packages=find_packages(),
    url='http://ghe-rss.roche.com/plsRED-Bioinformatics/',    
    author='Richard Dannebaum',
    author_email='Richard.Dannebaum@roche.com',
    description='A object for constructing string similarity networks, used for identifying groups of similar sequence(s) from sequencing datasets. Useful for clustering similar sequences or defining UMI families for PCR deduplication',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        'python-levenshtein',
        'pandas'
    ]
)
