// export PATH=/x86_64-w64-mingw32.static.posix/bin:$PATH
#include <iostream>
#include "nlohmann/json.hpp"
#include "hru.h"
#include <fstream>
using json = nlohmann::json;

int main() {
  // input files
  std::fstream File;
  File.open(R"(SwindaleModel.geojson)",std::ios::in);

  json Doc(json::parse(File));

  //std::cout << Doc.dump();
  //std::cout << "\n";
  
  auto nhru{ Doc["features"].size() };

  std::cout << Doc["features"].size();
  std::cout << "\n";

  std::vector<hru> hrus(nhru);

  for(auto ii=0; ii < nhru; ii++){
    //std::cout << Doc["features"].at(ii)["properties"]; //.get_to(hrus.at(ii));
    //std::cout << "\n";
    Doc["features"].at(ii)["properties"].get_to(hrus.at(ii));
  }

  // Players.at(0).Level[0] += 1;


  std::cout << hrus.at(100).uid;
  // std::cout << "\n";
  
  // json DocOut(Players);
  // std::cout << DocOut.dump();
  // std::cout << "\n";

  // for(int ii=0; ii < nPlayer; ii++){
  //   std::cout << "Name: " << Players.at(ii).Name;
  //   std::cout << "\nLevel: " << Players.at(ii).Level[0];
  // }

  // std::fstream outFile;
  // outFile.open(R"(output.json)", std::ios::out);
  // outFile << DocOut;
  // outFile.close();
}
