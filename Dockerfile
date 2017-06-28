FROM andrewosh/binder-base

ADD . $HOME/notebooks

USER root

RUN chown -R main:main $HOME/notebooks

USER main

WORKDIR $HOME/notebooks

RUN conda env create -n binder python=3.6

ENV PATH /home/main/anaconda/envs/binder/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR $HOME/notebooks/py_gmm

RUN /bin/bash -c "source activate binder && sh f2py.sh"

WORKDIR $HOME/notebooks

RUN /bin/bash -c "source activate binder && ipython kernelspec install-self --user"

iscarding /home/main/anaconda2/bin from PATH
prepending /home/main/anaconda2/envs/binder/bin to PATH
[TerminalIPythonApp] WARNING | Subcommand `ipython kernelspec` is deprecated and will be re
moved in future versions.
[TerminalIPythonApp] WARNING | You likely want to use `jupyter kernelspec` in the future
[InstallNativeKernelSpec] WARNING | `jupyter kernelspec install-self` is DEPRECATED as of 4
.0. You probably want `ipython kernel install` to install the IPython kernelspec.
[InstallNativeKernelSpec] Installed kernelspec python2 in /home/main/.local/share/jupyter/k
ernels/python2
