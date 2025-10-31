![DisPerSE](https://github.com/thierry-sousbie/DisPerSE/blob/master/manual/web/images/logo.png "DisPerSE")
Version 0.9.25  
Copyright(c) 2011- by Thierry Sousbie. All rights reserved.  
Author: Thierry Sousbie - tsousbie[a]gmail[dot]com  

DisPerSE stands for "Discrete Persistent Structures Extractor" and its main purpose is the automatic identification of persistent topological features such as peaks, voids, walls and in particular filamentary structures within sampled distributions in 2D, 3D, and possibly more ...

Although it was initially developed with cosmology in mind (for the study of the properties of filamentary structures in the so called comic web of galaxy distribution over large scales in the Universe), the present version is quite versatile and should be useful for any application where a robust structure identification is required, for segmentation or for studying the topology of sampled functions (like computing persistent Betti numbers for instance).

DisPerSE is able to deal directly with noisy datasets using the concept of persistence (a measure of the robustness of topological features) and can work indifferently on many kinds of cell complex (such as structured and unstructured grids, 2D manifolds embedded within a 3D space, discrete point samples using delaunay tesselation, Healpix tesselations of the sphere, ...). The only constraint is that the distribution must be defined over a manifold, possibly with boundaries. 

For more information, read these two scientific articles: [Theory](http://adsabs.harvard.edu/abs/2011MNRAS.414..350S "DisPerSE Theory article") and [Application](http://adsabs.harvard.edu/abs/2011MNRAS.414..384S "DisPerSE Application article").

You can also visit [DisPerSe](http://www2.iap.fr/users/sousbie/web/html/indexd41d.html) website at: http://www2.iap.fr/users/sousbie/web/html/indexd41d.html

---

# Instalation options

If a no-visualization setup is acceptable, a practical option could be use Docker.

The following docker image provides a full **no GUI** DisPerSE build (no `pdview` available)
``` bash
docker pull franucles/disperse:latest
```

It can be executed with
``` bash
docker run -it franucles/disperse:latest
```

If a complete instalation is required, the below instructions must be followed.

# Installing Required Packages

To avoid version conflicts with system packages, we will use a virtual environment.

Create a virtual environment for DisPerSE:

```bash
mamba env create -n disperse
```

The environment is now created, but we need to set a couple of environment variables for everything to work properly.

First, we need the path to the dynamic libraries. For example, run:

```bash
find /usr -name libmvec.so.1
```

Take the directory from the result (in my case `/usr/lib/x86_64-linux-gnu`) and add it with:

```bash
conda env config vars set LDFLAGS="-L/usr/lib/x86_64-linux-gnu" LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
```

> ℹ️ This uses `conda` instead of `mamba` because `mamba` doesn’t support `env config vars set`.

Now we can install the required and optional packages. I installed all of them to use every tool that DisPerSE offers.

---

## gcc

The most important tool since it is used to compile the code. Most systems already include gcc, but to be safe:

```bash
mamba install gcc
```

This will install the latest version. I got **gcc 15.1.0**, which works fine (minimum required is 4.3).

> ⚠️ Newer versions should work without issues. Tested successfully with gcc 15.1.0.

---

## CMake

Needed to generate Makefiles for building DisPerSE.

```bash
mamba install cmake
```

This installs the latest version.  
The repository requires version 2.8, but this is very outdated. I installed **CMake 4.1.1**.  

Since support for versions <3.5 was removed ([release log](https://cmake.org/cmake/help/latest/release/4.0.html#deprecated-and-removed-features)), some `CMakeLists` files must be modified. These changes are already included in [my repository](https://github.com/FranUcles/DisPerSE).

---

## GSL

A math library:

```bash
mamba install gsl
```

Installed version: **2.8** (no compatibility issues).

---

## CGAL

Optional. Used for Delaunay tessellation:

```bash
mamba install cgal
```

The repo requires 3.7, but I got **6.01**, which works fine.

---

## CFitsIO

Optional. For reading FITS files:

```bash
mamba install cfitsio
```

Installed version: **4.6.2**.

---

## SDL and SDL-image

Optional. Needed to insert images when using **pdview**:

```bash
mamba install sdl sdl_image
```

Installed version: **1.2.68**.

---

## Qt

Optional. Provides windows for **pdview** graphics:

```bash
mamba install qt
```

The repo requires Qt4, but I installed **Qt5**.  
This requires small changes in the `CMakeLists`, which are already included in [my repository](https://github.com/FranUcles/DisPerSE).

---

## MathGL

Optional, but required if you want to use **pdview**.  
It handles 2D/3D graphics and large datasets. Documentation: [MathGL official website](https://mathgl.sourceforge.net/).

Not available via conda or apt — must be built manually:

1. Download from the [official download page](https://mathgl.sourceforge.net/doc_en/Download.html).  
2. Extract and run:

```bash
cmake -Denable-qt5=ON -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .
```

- `-Denable-qt5=ON` → enables Qt5 support for DisPerSE  
- `-DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX` → installs inside the conda environment  

> ⚠️ Omit `-DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX` if you want a system-wide installation.  

3. Compile and install:

```bash
make -j$(nproc)
sudo make install -j$(nproc)
```

The official repo uses MathGL 1.11, but modern versions (e.g., **8.0.3**) are incompatible. I modified the DisPerSE code to work with the latest MathGL.

---

# Compiling and Installing DisPerSE

With all dependencies installed, we can compile the code:

1. Generate Makefile:

```bash
cmake .
```

2. Compile:

```bash
make -j$(nproc)
```

3. Install (binaries go into the `bin` folder):

```bash
make install -j$(nproc)
```

This installs between 2 and 5 executables:

```
* delaunay_2D (Optional)
* delaunay_3D (Optional)
* pdview (Optional)
* skelconv
* netconv
* mse
```
