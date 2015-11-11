FROM andrewosh/binder-base

USER main

ADD . $HOME/notebooks

USER root

RUN chown -R main:main $HOME/notebooks

USER main

RUN ls $HOME

RUN ls $HOME/notebooks
