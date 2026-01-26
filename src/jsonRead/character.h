#pragma once
#include "nlohmann/json.hpp"
#include <string>
using json = nlohmann::json;

class Character {
public:
  std::string Name;
  std::vector<int> Level;
};

void to_json(json& j, const Character& C) {
  j = json{{"Name", C.Name},
           {"Level", C.Level}};
}

void from_json(const json& j, Character& C) {
  j.at("Name").get_to(C.Name);
  j.at("Level").get_to(C.Level);
}
