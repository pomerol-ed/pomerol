#ifndef __INCLUDE__OUTPATH__
#define __INCLUDE__OUTPATH__

#include "config.h"
#include <string>

class output_handle
// A class to handle all output directory and file structure
{
  
  string path_str; // Output directory path

  void clean();

public:
  output_handle(){};
  output_handle(std::string Path);
  const string& path();
  string fullpath();
};

#endif // endif :: ___OUTPATH_H___
