#### tutorial to compile your code and link with pomerol
Here is a sample external compilation project to link to pomerol, consisting of 1 source file (example2site.cpp)

In order to compile it - run

```
cmake -Dpomerol_DIR=/where/pomerol/is/installed/share/pomerol /path/to/this/tutorial 
make
```
or alternatively add `/where/pomerol/is/installed/share/pomerol` to your `CMAKE_MODULE_PATH`. The latter is done automatically, when lmod is used with provided `pomerol.lmod`. 
