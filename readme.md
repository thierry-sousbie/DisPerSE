![DisPerSE](https://github.com/thierry-sousbie/DisPerSE/blob/master/manual/web/images/logo.png "DisPerSE")

DisPerSE stands for "Discrete Persistent Structures Extractor" and its main purpose is the automatic identification of persistent topological features such as peaks, voids, walls and in particular filamentary structures within sampled distributions in 2D, 3D, and possibly more ...

Although it was initially developed with cosmology in mind (for the study of the properties of filamentary structures in the so called comic web of galaxy distribution over large scales in the Universe), the present version is quite versatile and should be useful for any application where a robust structure identification is required, for segmentation or for studying the topology of sampled functions (like computing persistent Betti numbers for instance).

DisPerSE is able to deal directly with noisy datasets using the concept of persistence (a measure of the robustness of topological features) and can work indifferently on many kinds of cell complex (such as structured and unstructured grids, 2D manifolds embedded within a 3D space, discrete point samples using delaunay tesselation, Healpix tesselations of the sphere, ...). The only constraint is that the distribution must be defined over a manifold, possibly with boundaries. 

More information avilable at: http://www2.iap.fr/users/sousbie/web/html/indexd41d.html

REQUIREMENTS
------------

Any recent version of gcc should work, but gcc4.3.3+ is required in order to enable multithreading. The program was not tested with other compilers.

The following libraries/programs are used:

    Cmake 2.8     (Required for compilation)
    GSL	      (Required)
    CGAL 3.7      (Optional, for delaunay triangulation/DTFE)
    CFitsIO       (Optional, for reading FITS and Healpix maps)
    SDL/SDL-image (Optional, for loading jpg,bmp,...)
    mathGL 1.10+  (Optional, for visualizing persistence diagrams with pdview)
                  -> A compatible version can be found in the 'external/' subdir
    Qt4    	      (Optional, for visualizing persistence diagrams with pdview)

An appropriate version of mathGL can be found int the `external/` subdir.

COMPILATION
-----------

Suppose the source package is uncompressed in ${DISPERSE_SRC}, then go to the ${DISPERSE_SRC}/build directory and run cmake:
   
     cd ${DISPERSE_SRC}
     cd build
     cmake ../ 

This will check the configuration and generate the Makefile. Read the output to know which library were found, which were not, and how to specify their path (option -D{LIBNAME}_DIR="path/to/library/" where {LIBNAME} may be QT,GSL,SDL,MATHGL or CGAL).

Then just compile with:

     make
  or 
     make -j N (where N is the number of processors to use)


INSTALLATION
------------

Typing :

     make install

will install the programs in the "CMAKE_INSTALL_PREFIX" direction. By default, CMAKE_INSTALL_PREFIX=${DISPERSE_SRC}/bin, and the program is installed in the 'bin' subdirectory of the source package.
You can choose a different location with:

     cmake [...] -DCMAKE_INSTALL_PREFIX=PATH/TO/INSTALL/DIR

Two to five executables should be created in ${DISPERSE_SRC}/bin, depending on which libraries were found:

    * delaunay_2D (Optionnal)
    * delaunay_3D (Optionnal)
    * pdview (Optionnal)
    * skelconv
    * netconv
    * mse


TROUBLESHOOTING
---------------

LIBRAIRIES:

In case you need to install a library in a non standard path, you can specify its location like this:
   
     cmake [...] -D{LIBNAME}_DIR=path/to/library

where "LIBNAME" is the uppercase name of the library (ie QT,GSL,SDL,MATHGL or CGAL) and "path/to/library" points to the path where the library is installed (i.e. "path/to/library/lib" should contain the library itself and "path/to/library/includes" the header files).

COMPILE ERRORS:

Multithreading is enabled by default but so far, it works only with gcc4.3.3+. If you have an older version, the compiler will complain and return errors. You can choose to enable/disable multithreading with:
 
     cmake [...] -DUSE_THREADS=true/false
