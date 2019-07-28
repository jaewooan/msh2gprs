/************************************************************************/
/* DISCRETE FEATURE MODEL                                               */
/* by Mohammad Karimi-Fard (karimi@stanford.edu)                        */
/* March 2007, Stanford, CA.                                            */
/*                                                                      */
/* December 2012. Modifications by Timur Garipov                        */
/* Geomechanical interface. Convert Tetgen data to Karimi data           */
/************************************************************************/
#pragma once

#include <FlowData.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <algorithm>
#include <math.h>
#include <iterator>
#include <vector>
#include <set>
#include <time.h>


enum ConnectionType : int
{
  matrix_matrix = 1,
  matrix_fracture = 2,
  fracture_fracture = 3
};

namespace flow
{


class CalcTranses
{
public:
   CalcTranses();
  ~CalcTranses();
  void compute_flow_data();
  static void save_output(const FlowData    & data,
                          const std::string & output_dir);
  // this guy writes text output in a series of files
  void writeOutputFiles(const std::string & output_path) const;
  void extractData(FlowData & data) const;
  void init();


public:
  int NbNodes;
  int NbPolygons;
  int NbPolyhedra;
  int NbZones;
  int NbFracs;
  int NbOptions;
  // coordinates
  std::vector<double>	X,Y,Z;
  // faces
  std::vector<std::vector<std::size_t> > vvVFaces;
  std::vector<int> vCodePolygon;
  // elements
  std::vector<std::vector<std::size_t> > vvFNodes;
  std::vector<int> vCodePolyhedron;
  //properties
  std::vector<int> vZoneCode;
  std::vector<double> vZVolumeFactor;
  std::vector<double> vTimurConnectionFactor;
  std::vector<double> vZPorosity;
  std::vector<int> vZPermCode;
  std::vector<double> vZPermeability;
  std::vector<double> vZConduction;

protected:
  void ProjectionA(double mx,double my,double mz,
                   double px,double py,double pz,
                   double ix,double iy,double iz,
                   double nx,double ny,double nz,
                   double *hx,double *hy,double *hz);
  void ProjectionB(double mx,double my,double mz,
                   double ix,double iy,double iz,
                   double ux,double uy,double uz,
                   double *hx,double *hy,double *hz);

  void ComputeBasicGeometry();
  void ComputeControlVolumeList();
  void PrepareConnectionList();
  void ConstructConnectionList();
  void VolumeCorrection();
  void ComputeContinuityNode();
  void ComputeDirectionalPermeability();
  void ComputeTransmissibilityPart();
  void ComputeTransmissibilityList();

protected:
// SimData * pSim;

  std::size_t	NbCVs,NbVolumes,NbInterfaces,NbEquations,NbFeatures,NbCF,NbConnections,NbIntersections;
  std::size_t	NbTransmissibility,NbMetric;
  int		OptionVC,OptionGO,OptionMC;

  ///// Grid information /////
  std::vector<std::vector<std::size_t>>		FNodes;

  std::vector<int>		CodePolygon;
  std::vector<int>		CodePolyhedron;

  ///// Control volume numbers /////
  std::vector<int>		EQF;
  std::vector<int>		EQV;

  ///// Additional polygon information /////
  std::vector<double> face_area;
  // center of mass
  std::vector<double> polygon_center_x, polygon_center_y, polygon_center_z;
  std::vector<double> face_normal_x, face_normal_y, face_normal_z;		// normal unit vector

  ///// Additional polyhedron information /////
  // geometric cell volume
  std::vector<double> polyhedron_volume;
  // center of mass
  std::vector<double> polyhedron_center_x,
                      polyhedron_center_y,
                      polyhedron_center_z;

  ///// Definition of the control volumes /////
  double Tolerance;

  // control volumes data
  std::vector<int>		CVType;
  std::vector<int>		CVZone;
  std::vector<double>	CV_center_x, CV_center_y, CV_center_z;
  std::vector<double>		CVVolume;  // volume of control volume (including volume facetor)

  // zone data
  std::vector<int>		ZoneCode;
  std::vector<double>	ZVolumeFactor;
  std::vector<double>	ZPorosity;
  std::vector<int>		ZPermCode;
  std::vector<std::vector<double>>	ZPermeability;
  std::vector<std::vector<double>>	ZConduction;
  double		K1,K2,K3,K4,K5,K6;
  double		Kx,Ky,Kz;

///// Definition of the connections /////

  std::vector<int>		connection_type;
  std::vector<int>		n_connection_elements;
  // connection control volumes
  std::vector<std::vector<int>>		ConCV;
  std::vector<std::vector<double>>		ConTr;

  std::vector<std::vector<double>>		ConGeom;
  std::vector<std::vector<double>>		ConMult;

  std::vector<std::vector<double>>	ConArea;
  std::vector<std::vector<double>>	ConPerm;
  /// same as con_normal, idk why Mo needs this
  std::vector<double>		ConP1x, ConP1y, ConP1z;
  std::vector<double>		ConP2x, ConP2y, ConP2z;
  // centers of connections
  std::vector<double>		con_center_x, con_center_y, con_center_z;
  // connection normals
  std::vector<double>		con_normal_x, con_normal_y, con_normal_z;
  // projection of connection center onto line connecting two control volumes
  std::vector<double>		con_center_proj_x, con_center_proj_y, con_center_proj_z;

///// Transmissibility List /////

  std::vector<int>		iTr, jTr;
  std::vector<double>	Tij,SumTr, TConductionIJ;

  int		NbEdges;	// Actives
  std::vector<int>		ListV1, ListV2;
  std::vector<int>		ListE1, ListE2, ListF;

  double		TotalVolume,FaceArea,delta;
  double		LocalDistance,t,xt,yt,zt;

  clock_t		t1,t2,t3,t4,t5,t6,t7,t8,t9;
  clock_t		Deb_Computing,Fin_Computing;

 public:
  double fracporo;
  static constexpr double transmissibility_conversion_factor =
      0.0085267146719160104986876640419948;
};

}
