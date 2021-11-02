from distutils.core import setup

setup(
    name = 'maia',
    packages = ['maia'],
    version = '0.1',  # Ideally should be same as your GitHub release tag varsion
    description = 'Microwave passive filter design assistant tool',
    #long_description = readme,
    author = 'hirp7',
    author_email = '',
    url = 'https://github.com/hirp7/Maia',
    download_url = '',
    keywords = ['RF', 'Microwave Engineering','Filter Design'],
    classifiers = [],
    install_requires = [
        'numpy',
        'scipy',
        'scikit-rf'],
)