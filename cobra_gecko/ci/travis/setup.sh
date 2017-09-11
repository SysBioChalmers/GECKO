#!/usr/bin/env bash

set -e

dependencies=(pip setuptools wheel tox twine)

if [[ "${TRAVIS_OS_NAME}" == "osx" ]]; then
    brew update
    if [[ "${TOXENV}" == "py27" ]]; then
        brew upgrade python
        brew link --overwrite python
        pip2 install -U ${dependencies[@]}
    else
        brew install python3
        brew link --overwrite python3
        pip3 install -U ${dependencies[@]}
    fi
else
    pip install -U ${dependencies[@]}
fi
