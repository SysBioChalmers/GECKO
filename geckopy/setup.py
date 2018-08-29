#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
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
    name='geckopy',
    version='1.3.1',
    description="Methods for using the GECKO model with cobrapy",
    long_description=readme + '\n\n' + history,
    author="Benjamin Sanchez",
    author_email='bensan@chalmers.se',
    url='https://github.com/SysBioChalmers/GECKO/tree/master/geckopy',
    packages=[
        'geckopy',
    ],
    package_dir={'geckopy':
                 'geckopy'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT License",
    zip_safe=False,
    keywords='geckopy',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
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
