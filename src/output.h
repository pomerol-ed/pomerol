#ifndef __INCLUDE__OUTPATH__
#define __INCLUDE__OUTPATH__
#include "config.h"
#include <boost/filesystem.hpp>
#include <string>

namespace bf = boost :: filesystem;

class output_handle
// A class to handle all output directory and file structure
{
  
  string path_str; // Output directory path

  bf::path path_;

  void clean();

public:
  output_handle(){};
  output_handle(std::string Path);
  const string& path();
  string fullpath();
};

void progressbar(int percent);
#endif // endif :: ___OUTPATH_H___
