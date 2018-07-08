from setuptools import setup, find_packages
from configparser import ConfigParser

# Get some values from the setup.cfg
conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items('metadata'))


AUTHOR = metadata.get('author', '')
AUTHOR_EMAIL = metadata.get('author_email', '')
DESCRIPTION = metadata.get('description', '')
KEYWORDS = metadata.get('keywords', '')
LICENSE = metadata.get('unknown')
LONG_DESCRIPTION = metadata.get('long_description', '')
PACKAGENAME = metadata.get('package_name', 'packagename')


setup(name=PACKAGENAME,
      version='0.0.1',
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      packages=find_packages(),
      keywords=KEYWORDS,
      setup_requires=['pytest-runner'],
      classifiers=[
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      )
