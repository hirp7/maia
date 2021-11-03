#from distutils.core import setup
from setuptools import setup, find_packages
setup(
    name = 'maia',
    version = '0.1',  # Ideally should be same as your GitHub release tag varsion
    description = 'Microwave passive filter design assistant tool',
    #long_description = readme,
    author = 'hirp7',
    #author_email = '',
    packages = find_packages(exclude=('maia'))
    url = 'https://github.com/hirp7/Maia',
    download_url = '',
    keywords = ['RF', 'Microwave Engineering','Filter Design'],
    classifiers = [],
    license = license,
    install_requires = [
        'numpy',
        'scipy',
        'scikit-rf']
)
