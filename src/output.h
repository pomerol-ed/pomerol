
#ifndef __INCLUDE_OUTPAT_H
#define __INCLUDE_OUTPAT_H

#include "Misc.h"

class output_handle
// A class to handle all output directory and file structure
{
  
  std::string path_str; // Output directory path

  void clean();

public:
  output_handle(){};
  output_handle(std::string Path);
  const std::string& path();
  std::string fullpath();
};

#endif // endif :: __INCLUDE_OUTPAT_H
