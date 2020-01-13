#pragma once

#include "SimData.hpp"
#include "mesh/Mesh.hpp"
#include "VTKWriter.hpp"

namespace gprs_data
{

class OutputDataVTK
{
 public:
  OutputDataVTK(const SimData & sim_data, const VTKOutputConfig config);
  void write_output(const std::string & output_path);

 private:
  void save_reservoir_data(const std::string & fname);
  void save_dfm_data(const std::string & fname);
  void save_edfm_data(const std::string & fname);
  // size is different for mech and flow
  // void saveMultiScaleSupport(const multiscale::MultiScaleOutputData & ms,
  //                            const std::size_t                        size,
  //                            const std::string                      & prefix,
  //                            std::ofstream                          & out);

  const SimData & m_data;
  const mesh::Mesh & m_grid;
  VTKOutputConfig m_config;
};

}