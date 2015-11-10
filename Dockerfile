FROM andrewosh/binder-base
ADD gevero/py_gmm_binder4/environment.yml environment.yml
RUN conda env create -n binder
