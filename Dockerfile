FROM andrewosh/binder-base
ADD environment.yml environment.yml
RUN conda env create -n binder
