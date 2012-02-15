#include "Lattice.h"
#include <fstream>

namespace Pomerol{

std::ostream& operator<<(std::ostream& output, const Lattice::Site& out)
{
    output << "Site \"" << out.label << "\", " << out.OrbitalSize << " orbital" << ((out.OrbitalSize>1)?"s":"") << ", " << out.SpinSize << " spin" << ((out.SpinSize>1)?"s":"") << ".";
	return output;
}

Lattice::Lattice(){
};

JSONLattice::JSONLattice(){
root = new Json::Value;
};


int JSONLattice::readin(const std::string &filename)
{
  Json::Reader reader;
  std::ifstream in;
  in.open(filename.c_str());
  try
  {
    bool parsingSuccessful = reader.parse( in, *root );
    if ( !parsingSuccessful )
  	{
		std::cout  << "Failed to parse configuration\n";
		std::cout << reader.getFormatedErrorMessages();
        return 1;
  	}
  }
  catch (std::exception ErrorException)
  	{
		std::cout << ErrorException.what() << std::endl;
		exit(1);
  	}
  in.close();
  readSites((*root)["Sites"]);
  return 0;
};

void JSONLattice::readSites(Json::Value &JSONSites)
{
    Log.setDebugging(true);
    for (Json::Value::iterator it=JSONSites.begin(); it!=JSONSites.end(); ++it){
        DEBUG(*it);
        std::string label = it.key().asString();
        unsigned short OrbitalSize = (*it)["OrbitalSize"].asInt();
        DEBUG(OrbitalSize);
        DEBUG((*it)["U"].asDouble());
        }
    exit(0);
};

} // end of namespace Pomerol
