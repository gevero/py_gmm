FROM andrewosh/binder-base

MAINTAINER Andrew Osheroff <andrewosh@gmail.com>

# install conda environment
RUN /bin/bash -c "cd /home/giovi/notebooks"
RUN conda env create -n py3 -f environment.yml
RUN source activate py3

USER root

# Add gfortran
RUN apt-get update
RUN apt-get install -y gfortran

USER main

# Install py_gmm
ENV OPENBLAS_NUM_THREADS 1
RUN sh f2py.sh
