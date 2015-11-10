FROM andrewosh/binder-base
ADD repo/environment.yml environment.yml
RUN conda env create -n binder
