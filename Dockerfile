# 0) calling the base image
FROM gcr.io/generic-notebooks/binder-base

# 1) adding my environment.yml file
ADD repo/environment.yml environment.yml

# 2) replicating my environment
RUN conda env create -n binder

# 3) fixing the path
RUN echo "export PATH=/home/main/anaconda/envs/binder/bin/,$PATH" >> ~/.binder_start

# 4) installing jupyter (but maybe it is already there)
RUN conda install -n binder jupyter

# 5) activating environment and registrating kernels
RUN /bin/bash -c "source activate binder && ipython kernelspec install-self --user"

# 6) adding my notebooks home
ADD repo $HOME/notebooks

# 7) changing ownership of the folder
USER root
RUN chown -R main,main $HOME/notebooks
USER main

# 10) signing notebooks
RUN find $HOME/notebooks -name '*.ipynb' -exec ipython trust {} \;

# 11) defining working folder
WORKDIR $HOME/notebooks

# 12) compiling my code
sh f2py.sh
