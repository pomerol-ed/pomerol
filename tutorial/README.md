# Tutorial to compile your code and link with pomerol

Here is a sample external CMake project to link to pomerol.
It consists of just one source file, `example2site.cpp`.

In order to compile it, run

```shell
cmake -Dpomerol_DIR=<pomerol_installation_dir>/lib/cmake /path/to/this/tutorial
make
```

Visit https://pomerol-ed.github.io/pomerol/html/ for comprehensive documentation
of pomerol's API.
