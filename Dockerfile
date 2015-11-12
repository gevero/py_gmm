FROM andrewosh/binder-base

ADD . $HOME/notebooks

USER root

RUN chown -R main:main $HOME/notebooks

USER main

WORKDIR $HOME/notebooks

RUN conda env create -n binder

WORKDIR $HOME/notebooks/py_gmm

RUN /bin/bash -c "source activate binder && sh f2py.sh"

WORKDIR $HOME/notebooks

RUN echo \"export PATH=/home/main/anaconda/envs/binder/bin/:$PATH\" >> ~/.binder_start

RUN conda install -n binder jupyter

RUN /bin/bash -c "source activate binder && ipython kernelspec install-self --user"
