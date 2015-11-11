FROM andrewosh/binder-base

USER main

# ADD . $HOME/notebooks

# USER root

# RUN chown -R main:main $HOME/notebooks

USER main

ADD environment.yml environment.yml

RUN conda env create -n binder
