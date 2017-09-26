#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'cobra>=0.8.2'
]

test_requirements = [
     'pytest'
]

setup(
    name='cobra_gecko',
    version='0.0.1',
    description="Methods for using the GECKO model with cobrapy",
    long_description=readme + '\n\n' + history,
    author="Bejamin Sanchez",
    author_email='bensan@chalmers.se',
    url='https://github.com/Biosustain/cobra_gecko',
    packages=[
        'cobra_gecko',
    ],
    package_dir={'cobra_gecko':
                 'cobra_gecko'},
    include_package_data=True,
    install_requires=requirements,
    license="custom..",
    zip_safe=False,
    keywords='cobra_gecko',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
