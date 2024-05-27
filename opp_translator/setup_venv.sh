#!/bin/bash

#Do not run below command if virtualenv is already installed
python3 -m pip install --user virtualenv

if [[ -z $OPP_INSTALL_PATH ]]; then
    echo "Please set OPP_INSTALL_PATH before running this script"
    exit
fi

mkdir -p $OPP_INSTALL_PATH/../opp_translator/opp_venv

python3 -m venv $OPP_INSTALL_PATH/../opp_translator/opp_venv

source $OPP_INSTALL_PATH/../opp_translator/opp_venv/bin/activate

python3 -m pip install --upgrade pip

python3 -m pip install -r $OPP_INSTALL_PATH/../opp_translator/requirements.txt

python3 -m pip install --force-reinstall libclang==16.0.6

