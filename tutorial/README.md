# Tutorial to compile your code and link with pomerol

Here is a sample external CMake project to link to pomerol.
It consists of just one source file, `example2site.cpp`.

In order to compile it, run

```shell
cmake -Dlibcommute_DIR=<libcommute_installation_dir>/lib/cmake/libcommute \
      -Dpomerol_DIR=<pomerol_installation_dir>/share/pomerol \
      /path/to/this/tutorial
make
```

Visit https://pomerol-ed.github.io/pomerol/html/ for comprehensive documentation
of pomerol's API.
