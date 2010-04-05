#include "output.h"
#include "config.h"

#include <boost/version.hpp>

#include <iostream>
#include <sstream>

output_handle::output_handle (std::string Path)
{
  path_str = Path;

  path_ = bf::path (path_str);

  if (bf::exists(path_))
    {
 //      std::cout<<"Output directory\""<<path_ << "\" exists. Using it" << endl;
//       (*this).clean();
    }
  else 
    {
      std::cout<<"Creating directory \""<<path_ << "\" for output" << endl;
      bf::create_directory(path_);
    }
};

void output_handle::clean()
{  
  if (bf::exists(path_))
    {
       std::cout<<"Cleaning output directory\""<<path_ << "\"" << endl;
       bf::directory_iterator dir_iter(path_), dir_end;
       #if BOOST_VERSION >= 103600
         for(;dir_iter != dir_end; ++dir_iter) if (bf::is_regular_file(dir_iter->status())) bf::remove(*dir_iter);
       #else
         for(;dir_iter != dir_end; ++dir_iter) if (bf::is_regular(dir_iter->status())) bf::remove(*dir_iter);
       #endif
    }

}

const string& output_handle::path()
{
  return path_str;
}

string output_handle::fullpath()
{
  return (bf::system_complete(path_)).directory_string();
}

void progressbar(int percent)
/* Prints a progressbar */
{
  static std::stringstream bars;
  static int x = 0;
  string slash[4];
  slash[0] = "\\";
  slash[1] = "-";
  slash[2] = "/";
  slash[3] = "|";
//  bars << "|";
  //cout << "\r"; // carriage return back to beginning of line
//  cout << bars.str() << " " << slash[x] << " " << percent << " %"; // print the bars and percentage
  cout << percent << " " << flush;
  x++; // increment to make the slash appear to rotate
  if(x == 4)
  x = 0; // reset slash animation
}
