#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/container/map/map_fwd.hpp>
#include <boost/fusion/include/map_fwd.hpp>
#include <iostream>
#include <map>
#include <string>
#include "Config.h"

void Config::load(const std::string &filename)
{
  using boost::property_tree::ptree;
  ptree pt;

  read_xml(filename, pt);
  this->databaseInfo["number_of_errors"] = pt.get<std::string>("bch.number_od_errors");
  this->databaseInfo["polynomial_degree"] = pt.get<std::string>("bch.polynomial_degree");
  this->databaseInfo["correction_capabilities"] = pt.get<std::string>("bch.correction_capabilities");
  
}

std::map<std::string, std::string> Config::getDatabaseInfo(){
  return this->databaseInfo;
}
