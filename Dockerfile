FROM andrewosh/binder-base
ADD py_gmm_binder4/environment.yml environment.yml
RUN conda env create -n binder
