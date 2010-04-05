#include <iostream>
#include <fstream>
#include <json/json.h>
#include <complex>

using namespace std;

int main()
{
Json::Value root;   // will contains the root value after parsing.
Json::Reader reader;
std::ifstream in("example.json");
bool parsingSuccessful = reader.parse( in, root );
if ( !parsingSuccessful )
{
    // report to the user the failure and their locations in the document.
	std::cout  << "Failed to parse configuration\n";
	std::cout << reader.getFormatedErrorMessages();
    return 1;
}

// Get the value of the member of root named 'encoding', return 'UTF-8' if there is no
// such member.
std::string encoding = root.get("encoding", "UHTF-8" ).asString();
cout << encoding << endl;
// Get the value of the member of root named 'encoding', return a 'null' value if
// there is no such member.
int Nbit = root["N_bit"].asInt()
const Json::Value plugins = root["plug-ins"];
for ( int index = 0; index < plugins.size(); ++index )  // Iterates over the sequence elements.
   std::cout<< plugins[index].asString();


 double shit = root["JJ"].asDouble();
 cout << endl << shit << endl;

   
return 1;
}
