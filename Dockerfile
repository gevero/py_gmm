FROM jupyter/datascience-notebook

ADD . $HOME/notebooks

USER root

RUN apt-get install libopenblas-base

RUN chown -R main:main $HOME/notebooks

USER main

WORKDIR $HOME/notebooks

RUN conda env create -n binder python=3.6

ENV PATH /home/main/anaconda2/envs/binder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR $HOME/notebooks/py_gmm

RUN /bin/bash -c "source activate binder && sh f2py.sh"

WORKDIR $HOME/notebooks

# RUN /bin/bash -c "source activate binder && ipython kernelspec install-self --user"
