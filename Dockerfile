FROM andrewosh/binder-base

USER root

RUN ls .

RUN ls $HOME

USER main

RUN ls .

RUN ls $HOME

ADD . $HOME/notebooks

USER root

RUN chown -R main:main $HOME/notebooks

USER main

RUN ls $HOME

RUN ls $HOME/notebooks
