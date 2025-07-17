#pragma once

#include <string>
#include <vector>

#include "dimensions.hpp"
#include "instance.hpp"

std::vector<Task> randomOrder(const sualbsp::Instance &);

unsigned short assign(const sualbsp::Instance &, const std::vector<Task> &, std::vector<short unsigned> &);
unsigned short assign(const sualbsp::Instance &, const std::vector<Task> &);

unsigned short sample(sualbsp::Instance &, unsigned short);

unsigned short assignRule(sualbsp::Instance &, unsigned short, std::string);
unsigned short assignRuleOpt(sualbsp::Instance &, unsigned short, std::string);
unsigned short assignRuleAlt(sualbsp::Instance &, sualbsp::Instance &, unsigned short, std::string);
unsigned short assignRuleOptAlt(sualbsp::Instance &, sualbsp::Instance &, unsigned short, std::string);