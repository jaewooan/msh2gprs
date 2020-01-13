#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "discretization/DoFNumbering.hpp"
#include "muparser/muParser.h" // parser for user-defined expressions for reservoir data

namespace gprs_data {

class CellPropertyManager
{
 public:
  CellPropertyManager(const CellPropertyConfig & cell_properties,
                      const std::vector<DomainConfig> & domain_configs,
                      SimData & data);
  void generate_properties();
  void map_mechanics_to_control_volumes(const discretization::DoFNumbering & dofs);

 private:
  void print_setup_message_();
  void assign_expressions_(const DomainConfig& domain,
                           std::vector<mu::Parser> & parsers,
                           std::vector<double> & vars);
  void evaluate_expressions_(const DomainConfig& domain,
                           std::vector<mu::Parser> & parsers,
                           std::vector<double> & vars);
  void build_permeability_function_();
  void build_porosity_function_();
  void build_flow_output_property_keys_();

  const CellPropertyConfig        & config;
  const std::vector<DomainConfig> & domains;
  SimData & m_data;
  // number of default variable in config
  // these variables are not output, so I don't save them
  // should be 3=x+y+z
  const size_t m_shift;
};

}  // end namespace gprs_data