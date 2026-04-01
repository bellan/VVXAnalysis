#!/bin/bash

# Please setup python 3 and ROOT into your environment first


init(){
    SUBSYSTEM=$1
    PACKAGE=$2

    CWD=$PWD
    if [ ${BASH_SOURCE[0]:0:1} == "/" ]; then
	FULLPATH=${BASH_SOURCE[0]}
    else
	FULLPATH=$PWD/${BASH_SOURCE[0]}
    fi
    cd ${FULLPATH/%env_standalone.sh/}/..
    echo $PWD
    cd ${CWD}/$SUBSYSTEM/$PACKAGE
    echo $PWD
    
    if [ ! -d build ]; then
	if [ x${3} = 'xbuild' ]; then
	    mkdir -p build/lib/python/$SUBSYSTEM
	    ln -s ../../../../python build/lib/python/$SUBSYSTEM/$PACKAGE
	    echo "Build directory created, please source again standalone/env_standalone.sh without the build argument."
	else
	    echo "Build directory is not yet present, please source again standalone/env_standalone.sh with the build argument."
	fi
    else
	if [ x${3} = 'xbuild' ]; then
	    echo "Build directory is already present, please source again standalone/env_standalone.sh without the build argument."
	else
	    find build/lib/python python -type d -execdir touch '{}/__init__.py' \;
	    export NANOAODTOOLS_BASE=${PWD}
	    export PYTHONPATH=${NANOAODTOOLS_BASE}/build/lib/python:${PYTHONPATH}
	    echo $1/$2 ": Standalone environment set."
	fi
    fi
    cd $CWD
}

