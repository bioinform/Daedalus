from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='daedalus_db',
    version='0.0.1',
    packages=["daedalus_db"],
    url='http://ghe-rss.roche.com/plsRED-Bioinformatics/Daedalus',
    author='Richard Dannebaum',
    author_email='richard.dannebaum@roche.com',
    description='A database to store sequencing run information',
    install_requires=[
        'pandas',
        'sqlalchemy'
    ]
)

