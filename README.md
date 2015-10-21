# py_matrix


A Fortran 90 implementation of the **Generalized Multiparticle Mie** method. If you find py_gmm useful for generating results included in publications, please consider citing the following paper:

Pellegrini G., Mattei G., Bello V. and Mazzoldi P. (2007) **Interacting metal nanoparticles: Optical properties from nanoparticle dimers to core-satellite systems** *Materials Science and Engineering: C 27:1347–1350*. [doi:10.1016/j.msec.2006.07.025](http://dx.doi.org/10.1016/j.msec.2006.07.025)

## Installation

Installing **py_gmm** should be fairly easy. First of all you need a python distribution on your system. The simplest thing to do, for any operating system (**Windows**, **OSX** or **Linux**), would be to install [Anaconda](https://store.continuum.io/cshop/anaconda/), a beautiful, free, and easy to use python distribution: the relevant installations instructions are found [here](http://docs.continuum.io/anaconda/install.html). Then you must also install The [Jupyter Notebook](http://ipython.org/notebook.html), which is already bundled in Anaconda.

The second step would be installing a Fortran compiler. If you are not new to compilers, you're already covered here. If don't know much about compilers yet, my advice is to go for [gfortran](https://gcc.gnu.org/wiki/GFortran). While I strongly advise you to work under **Linux** (this code was developed under **Linux**), I will try to provide minimal informations for the installations under **WIndows** and **OSX**. Keep in mind that I never installed or tried my software under these platforms, so I offer no guarantee whatsoever. Nevertheless if you encounter any problem please open an **Issue** here on github and we'll try to work it out. Regarding installation on different platforms it goes as follows:

* **Windows:** Get [MinGW](http://www.mingw.org/wiki/howto_install_the_mingw_gcc_compiler_suite) and follow the instructions to install **gfortran** on your system;
* **OSX:** Get the right [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS) and if needed read the instructions [here](https://gcc.gnu.org/wiki/GFortranBinariesMacOS);
* **Linux:** Install **gfortran** with the package manager of your own distro: it should be no trouble at all.

The third step would be to get **py_gmm** itself. If you are familiar with [git](http://git-scm.com/) and [github](https://github.com/) you can simply clone the repository, otherwise just [download](https://github.com/gevero/py_gmm/archive/master.zip) the zipped version of the repo and unpack it wherever you like.

Then we finally come to compiling: if you are under **Linux** or **OSX** `sh f2py.sh` should be sufficient to build everything. Under **Windows** try and change `f2py.sh` to `f2py.bat` and see happens. If it does not work open an issue and we'll work it out.

## Usage

The best thing to do for the moment is to start from the **.ipynb** files in the [examples](https://github.com/gevero/py_gmm/tree/master/examples) folder. You can load them in your local **Jupyter Notebook** instance: they should give you a fair idea about how to proceed for your calculations. Each function in the code features a detailed documentation easily accessible with the `Shift-Tab` tool-tip shortcut from the **Jupyter Notebook** interface.
