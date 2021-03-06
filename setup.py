#from distutils.core import setup
import setuptools
import os, sys
from setuptools import setup, find_packages

setup(
    name = 'maia',
    version = '0.0.1',  # Ideally should be same as your GitHub release tag varsion
    description = 'Microwave passive filter design assistant tool',
    #long_description = readme,
    author = 'hirp7',
    author_email = 'izunyan@gmail.com',
    #packages = ['maia'],
    packages=find_packages(),
    url = 'https://github.com/hirp7/maia',
    #download_url = '',
    keywords = ['RF', 'Microwave Engineering','Filter Design'],
    #classifiers = [],
    license = 'MIT',
    install_requires = [
        'numpy',
        'scipy',
        'scikit-rf',
        'wheel',
        'schemdraw']
)


"""
#with open("README.md", "r", encoding="utf-8") as fh:
#    long_description = fh.read()

setuptools.setup(
    name='test',
    version='0.0.1',
    author='hirp7',
    author_email='izunyan7@gmail.com',
    description='Testing installation of Package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/hirp7/test',
    #project_urls = {
    #    "Bug Tracker": "https://github.com/mike-huls/toolbox/issues"
    #},
    license='MIT',
    packages=['test'],
    install_requires=['requests'],
)
"""
