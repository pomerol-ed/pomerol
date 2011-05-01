#include "output.h"

#include <sstream>

#include <sys/stat.h>
#include <sys/types.h>
#include <limits.h>
#include <stdlib.h>

output_handle::output_handle (std::string Path)
{
    path_str = Path;
    struct stat OutputDir;

    if (!stat(path_str.c_str(),&OutputDir) && S_ISDIR(OutputDir.st_mode))
    {
        std::cout<<"Output directory\""<< path_str << "\" exists. Using it" << std::endl;
    }
    else 
    {
        std::cout<<"Creating directory \""<< path_str << "\" for output" << std::endl;
        mkdir(path_str.c_str(),0777);
    }
};

const std::string& output_handle::path()
{
    return path_str;
}

std::string output_handle::fullpath()
{
    char RealPathBuf[PATH_MAX];
    if(realpath(path_str.c_str(),RealPathBuf))
        return std::string(RealPathBuf);
    else
        return "";
}
