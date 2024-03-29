FROM jupyter/datascience-notebook:ubuntu-20.04

ADD . $HOME/notebooks

USER root

RUN apt update


RUN apt -y install liblapack-dev liblapack3 libopenblas-base libopenblas-dev gfortran
RUN ln -s /opt/OpenBLAS/lib/libopenblas.so /usr/lib/libopenblas.so

RUN chown -R $NB_USER:users $HOME/notebooks

USER $NB_USER

# WORKDIR $HOME/notebooks

# RUN conda env create -n binder python=3.6
# RUN pip install plotly

# ENV PATH /home/main/anaconda2/envs/binder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR $HOME/notebooks/py_gmm

# RUN /bin/bash -c "source activate binder && sh f2py.sh"
RUN /bin/bash -c "sh f2py.sh"

WORKDIR $HOME/notebooks

# RUN /bin/bash -c "source activate binder && ipython kernelspec install-self --user"
