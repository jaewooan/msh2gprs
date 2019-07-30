#define _CRT_SECURE_NO_DEPRECATE
#include "transes.hpp"
#include "simdata.hpp"
#include <random>

namespace flow
{


CalcTranses::CalcTranses()
{}


void CalcTranses::init()
{
  //coordinates
  X.resize(NbNodes);
  Y.resize(NbNodes);
  Z.resize(NbNodes);

  //faces
  vvVFaces.resize(NbPolyhedra);
  //elements
  // vNbFNodes.resize(NbPolygons, 3);
  vvFNodes.resize(NbPolygons);
  vCodePolyhedron.resize(NbPolyhedra);
  //properties
  vZoneCode.resize(NbZones);
  vZVolumeFactor.resize(NbZones);
  vTimurConnectionFactor.resize(NbZones);
  vZPorosity.resize(NbZones);
  vZPermCode.resize(NbZones);
  vZPermeability.resize(NbZones * 3, 0.0);
  vZConduction.resize(NbZones * 3, 0.0);

  OptionVC = 0;
  OptionGO = 0;
  OptionMC = 0;
  Tolerance = 0.05;

  // vNbFNodes.resize(NbPolygons);
  vvFNodes.resize(NbPolygons);
  vCodePolygon.resize(NbPolygons);

  ListV1.resize(NbPolygons);
  ListV2.resize(NbPolygons);
}

CalcTranses::~CalcTranses()
{
  // for (std::size_t)
  // free(vNbFNodes);
  // // **FNodes -- need to remove subelements but i'm lazy
  // for ()
  // free(FNodes);
  // free(NbVFaces);
  // free (VFaces);
}
/********************************************************************/
void CalcTranses::
ProjectionA( double CVx,double CVy,double CVz,  // CV center
             double px,double py,double pz,     // CV1 - CV2
             double face_center_x,double face_center_y,double face_center_z, // face center
             double nx,double ny,double nz,    // face normal
             double *hx,double *hy,double *hz)
{
  // (face_center - cell_center) · face_normal /
  // (cell1_center - cell2_center) · face_normal
  const double t = ((face_center_x-CVx)*nx +
                    (face_center_y-CVy)*ny +
                    (face_center_z-CVz)*nz) / (px*nx + py*ny + pz*nz);
  *hx = CVx + t*px;
  *hy = CVy + t*py;
  *hz = CVz + t*pz;
}
/********************************************************************/
void CalcTranses::ProjectionB( double mx,  double my,  double mz,
                               double ix,  double iy,  double iz,
                               double ux,  double uy,  double uz,
                               double *hx, double *hy, double *hz)
{
    double  t;

    t = (ux*(mx-ix) + uy*(my-iy) + uz*(mz-iz)) / (ux*ux + uy*uy + uz*uz);

    *hx = ix + t*ux;
    *hy = iy + t*uy;
    *hz = iz + t*uz;
}
/********************************************************************/
void CalcTranses::ComputeBasicGeometry()
{
    // int i, j, k;
    // int k;
    double  xi, yi, zi, areatmp, volumetmp;
    double  ux, uy, uz,
        vx, vy, vz,
        nx,ny,nz,
        nl;

//////////////////////////////////////////////////////////////////
///// Polygon Area, Center of Mass, and Normal (Unit vector) /////
//////////////////////////////////////////////////////////////////

    face_area.resize(NbPolygons);
    polygon_center_x.resize(NbPolygons);
    polygon_center_y.resize(NbPolygons);
    polygon_center_z.resize(NbPolygons);
    face_normal_x.resize(NbPolygons);
    face_normal_y.resize(NbPolygons);
    face_normal_z.resize(NbPolygons);

    for (std::size_t i=0; i<NbPolygons; i++)
    {
      face_area[i] = polygon_center_x[i] = polygon_center_y[i] = polygon_center_z[i] = 0;
      face_normal_x[i] = face_normal_y[i] = face_normal_z[i] = 0;
      for (std::size_t j=1; j<vvFNodes[i].size()-1; j++)
      {
        /* compute face normal and offset w from first 3 vertices */
        ux = X[FNodes[i][j]] - X[FNodes[i][0]];
        uy = Y[FNodes[i][j]] - Y[FNodes[i][0]];
        uz = Z[FNodes[i][j]] - Z[FNodes[i][0]];
        vx = X[FNodes[i][j+1]] - X[FNodes[i][0]];
        vy = Y[FNodes[i][j+1]] - Y[FNodes[i][0]];
        vz = Z[FNodes[i][j+1]] - Z[FNodes[i][0]];
        // not normalized normal vector components
        nx = (uy*vz - vy*uz);
        ny = (vx*uz - ux*vz);
        nz = (ux*vy - uy*vx);
        areatmp = .5*sqrt(nx*nx + ny*ny + nz*nz);

        face_area[i] += areatmp;
        polygon_center_x[i] += areatmp*(X[FNodes[i][0]] + X[FNodes[i][j]] + X[FNodes[i][j+1]])/3.;
        polygon_center_y[i] += areatmp*(Y[FNodes[i][0]] + Y[FNodes[i][j]] + Y[FNodes[i][j+1]])/3.;
        polygon_center_z[i] += areatmp*(Z[FNodes[i][0]] + Z[FNodes[i][j]] + Z[FNodes[i][j+1]])/3.;

        face_normal_x[i] +=.5*nx;
        face_normal_y[i] +=.5*ny;
        face_normal_z[i] +=.5*nz;
      }

      polygon_center_x[i] = polygon_center_x[i] / face_area[i];
      polygon_center_y[i] = polygon_center_y[i] / face_area[i];
      polygon_center_z[i] = polygon_center_z[i] / face_area[i];

      face_normal_x[i] = face_normal_x[i] / face_area[i];
      face_normal_y[i] = face_normal_y[i] / face_area[i];
      face_normal_z[i] = face_normal_z[i] / face_area[i];

      nl = sqrt(face_normal_x[i]*face_normal_x[i] + face_normal_y[i]*face_normal_y[i] + face_normal_z[i]*face_normal_z[i]);

      face_normal_x[i] = face_normal_x[i] / nl;
      face_normal_y[i] = face_normal_y[i] / nl;
      face_normal_z[i] = face_normal_z[i] / nl;
    }

////////////////////////////////////////////////
///// Polyhedron Volume and Center of Mass /////
////////////////////////////////////////////////

    polyhedron_volume.resize(NbPolyhedra);
    polyhedron_center_x.resize(NbPolyhedra);
    polyhedron_center_y.resize(NbPolyhedra);
    polyhedron_center_z.resize(NbPolyhedra);

    for (std::size_t i=0; i<NbPolyhedra; i++)
    {
        polyhedron_volume[i] = polyhedron_center_x[i] = polyhedron_center_y[i] = polyhedron_center_z[i] = 0;
        xi = yi = zi = 0;

        // Defining a node inside the polyhedron
        // it's the centroid of faces
        for (std::size_t j=0; j<vvVFaces[i].size(); j++)
        {
          xi += polygon_center_x[vvVFaces[i][j]];
          yi += polygon_center_y[vvVFaces[i][j]];
          zi += polygon_center_z[vvVFaces[i][j]];
        }

        const double n_element_faces = vvVFaces[i].size();
        xi = xi / n_element_faces;
        yi = yi / n_element_faces;
        zi = zi / n_element_faces;

        for (std::size_t j=0; j<vvVFaces[i].size(); j++)
        {
          const std::size_t face = vvVFaces[i][j];
          const double h = face_normal_x[face] * (polygon_center_x[face]-xi) +
                           face_normal_y[face] * (polygon_center_y[face]-yi) +
                           face_normal_z[face] * (polygon_center_z[face]-zi);

          volumetmp = fabs(h*face_area[face]) / 3.;

          if (std::isnan(volumetmp))
          {
            std::cout << "found nan volume" << std::endl;
            std::cout << "Polyhedron = " << i << std::endl;
            std::cout << "face = " << face << std::endl;
            abort();
          }

          polyhedron_center_x[i] += (polygon_center_x[face] +
                                     .25*(xi-polygon_center_x[face])) * volumetmp;
          polyhedron_center_y[i] += (polygon_center_y[face] +
                                     .25*(yi-polygon_center_y[face])) * volumetmp;
          polyhedron_center_z[i] += (polygon_center_z[face] +
                                     .25*(zi-polygon_center_z[face])) * volumetmp;
          polyhedron_volume[i] += volumetmp;
        }

        polyhedron_center_x[i] = polyhedron_center_x[i] / polyhedron_volume[i];
        polyhedron_center_y[i] = polyhedron_center_y[i] / polyhedron_volume[i];
        polyhedron_center_z[i] = polyhedron_center_z[i] / polyhedron_volume[i];

    }
}
/********************************************************************/
void CalcTranses::ComputeControlVolumeList()
{
  // int j;

  CVType.resize(NbCVs);
  CVZone.resize(NbCVs);
  CVVolume.resize(NbCVs);
  CV_center_x.resize(NbCVs);
  CV_center_y.resize(NbCVs);
  CV_center_z.resize(NbCVs);

  for (std::size_t i=0; i<NbPolygons; i++)
    if (EQF[i] != -1) // Active polygon
    {
      std::size_t j = EQF[i];
      CVZone[j] = CodePolygon[i];
      CVType[j] = 2;  // Feature
      CV_center_x[j] = polygon_center_x[i];
      CV_center_y[j] = polygon_center_y[i];
      CV_center_z[j] = polygon_center_z[i];
      CVVolume[j] = face_area[i] * ZVolumeFactor[CVZone[j]];
    }

  for (int i=0;i<NbPolyhedra;i++)
  {
    if (EQV[i] != -1) // Active polyhedron
    {
      std::size_t j = EQV[i];
      CVZone[j] = CodePolyhedron[i];
      CVType[j] = 1;  // Volume
      CV_center_x[j] = polyhedron_center_x[i];
      CV_center_y[j] = polyhedron_center_y[i];
      CV_center_z[j] = polyhedron_center_z[i];
      CVVolume[j] = polyhedron_volume[i] * ZVolumeFactor[CVZone[j]];
    }
  }
}
/********************************************************************/
void CalcTranses::PrepareConnectionList()
{
    // int i,j,k,iswap;
    int iswap;

    ///// Direct construction of M-M and M-F connections /////
    for (std::size_t i=0; i<NbPolygons; i++)
    {
      ListV1[i] = -1;
      ListV2[i] = -1;
    }

    for (std::size_t i=0; i<NbPolyhedra; i++)
    {
      for (std::size_t j=0; j<vvVFaces[i].size(); j++)
      {
        if (ListV1[vvVFaces[i][j]] == -1)
          ListV1[vvVFaces[i][j]] = i;
        else
          ListV2[vvVFaces[i][j]] = i;
      }
    }

    ///// Construction of F-F needs sorting... less efficient /////

    ListE1.resize(NbEdges);
    ListE2.resize(NbEdges);
    ListF.resize(NbEdges);

    std::size_t k = 0;
    for (std::size_t i=0; i<NbPolygons; i++)
      if (CodePolygon[i] >= 0)
      {
        for (std::size_t j=0; j<vvFNodes[i].size()-1; j++)
        {
          if (FNodes[i][j] < FNodes[i][j+1])
          {
            ListE1[k] = FNodes[i][j];
            ListE2[k] = FNodes[i][j+1];
            ListF [k] = i;
          }
          else
          {
            ListE1[k] = FNodes[i][j+1];
            ListE2[k] = FNodes[i][j];
            ListF [k] = i;
          }
          k++;
        }

        if (FNodes[i][0] < FNodes[i][vvFNodes[i].size()-1])
        {
          ListE1[k] = FNodes[i][0];
          ListE2[k] = FNodes[i][vvFNodes[i].size()-1];
          ListF [k] = i;
        }
        else
        {
          ListE1[k] = FNodes[i][vvFNodes[i].size()-1];
          ListE2[k] = FNodes[i][0];
          ListF [k] = i;
        }

        k++;
      }

    ///// Two step sorting /////

    // First step

    for (std::size_t i=0; i<NbEdges-1; i++)
    {
      if (NbEdges == 0)
        break;

      for (std::size_t j=i+1; j<NbEdges; j++)
      {
        if (ListE1[i] > ListE1[j])
        {
          iswap     = ListE1[i];
          ListE1[i] = ListE1[j];
          ListE1[j] = iswap;

          iswap     = ListE2[i];
          ListE2[i] = ListE2[j];
          ListE2[j] = iswap;

          iswap     = ListF[i];
          ListF[i]  = ListF[j];
          ListF[j]  = iswap;
        }
      }
    }

    // Second step

    for (std::size_t i=0; i<NbEdges-1; i++)
    {
      if (NbEdges == 0)
        break;

      for (std::size_t j=i+1; j<NbEdges; j++)
      {
        if (ListE1[i] == ListE1[j] && ListE2[i] > ListE2[j])
        {
          iswap     = ListE2[i];
          ListE2[i] = ListE2[j];
          ListE2[j] = iswap;

          iswap    = ListF[i];
          ListF[i] = ListF[j];
          ListF[j] = iswap;
        }
      }
    }

    ///// Evaluation of the number of connections /////

    // M-M and M-F

    NbConnections = 0;
    for (std::size_t i=0; i<NbPolygons; i++)
    {
      // if (ListV1[i] < 0)
      if (CodePolygon[i] < 0)
      {
        if (ListV2[i] != -1)
          NbConnections++;
      }
      else
      {
        if (ListV2[i] == -1)
          NbConnections++;
        else
          NbConnections += 2;
      }

    }

    // F-F

    int i = 0;
    while (i < NbEdges-1)
    {

      std::size_t j = i + 1;
      while ((j < NbEdges) && (ListE1[i] == ListE1[j]) && (ListE2[i] == ListE2[j]))
        j++;

      if ((j-i) >= 2)
        NbConnections++;

      i = j;
    }
}
/********************************************************************/
void CalcTranses::ConstructConnectionList()
{
    int j,k;

    // std::cout << "NbConnections = "<< NbConnections << std::endl;

    connection_type.resize(NbConnections);
    n_connection_elements.resize(NbConnections);
    ConCV.resize(NbConnections);
    ConTr.resize(NbConnections);
    ConArea.resize(NbConnections);

    ConGeom.resize(NbConnections);
    ConMult.resize(NbConnections);

    ConPerm.resize(NbConnections);
    ConP1x.resize(NbConnections);
    ConP1y.resize(NbConnections);
    ConP1z.resize(NbConnections);

    ConP2x.resize(NbConnections);
    ConP2y.resize(NbConnections);
    ConP2z.resize(NbConnections);

    con_center_x.resize(NbConnections);
    con_center_y.resize(NbConnections);
    con_center_z.resize(NbConnections);

    con_normal_x.resize(NbConnections);
    con_normal_y.resize(NbConnections);
    con_normal_z.resize(NbConnections);


    NbTransmissibility = 0;

    k=0;  // connection index
    if (NbPolyhedra > 0)
    for (std::size_t i=0; i<NbPolygons; i++)
    {
      if (CodePolygon[i] < 0 && ListV2[i] >= 0)  // M-M
      {
        // std::cout << "M-M " << i << std::endl;
        connection_type[k] = ConnectionType::matrix_matrix;
        n_connection_elements[k] = 2;

        ConCV[k].resize(2);
        ConCV[k][0] = EQV[ListV1[i]];
        ConCV[k][1] = EQV[ListV2[i]];

        ConArea[k].resize(2);
        ConArea[k][0] = ConArea[k][1] = face_area[i];

        ConP1x[k] = face_normal_x[i];
        ConP1y[k] = face_normal_y[i];
        ConP1z[k] = face_normal_z[i];

        ConP2x[k] = face_normal_x[i];
        ConP2y[k] = face_normal_y[i];
        ConP2z[k] = face_normal_z[i];

        con_center_x[k] = polygon_center_x[i];
        con_center_y[k] = polygon_center_y[i];
        con_center_z[k] = polygon_center_z[i];

        con_normal_x[k] = face_normal_x[i];
        con_normal_y[k] = face_normal_y[i];
        con_normal_z[k] = face_normal_z[i];

        ConTr[k].resize(2);
        ConPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;
        NbTransmissibility++;
      }

      // if (ListV1[i] >= 0 && ListV2[i] < 0) // M-F
      if (CodePolygon[i] >= 0 && ListV2[i] < 0) // M-F
      {
        // std::cout << "M-F " << i << std::endl;
        connection_type[k] = ConnectionType::matrix_fracture;
        n_connection_elements[k] = 2;

        ConCV[k].resize(2);

        if (ListV1[i] < 0)
        {
          std::cout << "polygon " << i << "\t";
          std::cout << "wrong ListV2 index " << ListV1[i] << std::endl;
          std::cout << "aborting" << std::endl;
          exit(0);
        }

        ConCV[k][0] = EQV[ListV1[i]];
        ConCV[k][1] = EQF[i];

        ConArea[k].resize(2);
        ConArea[k][0] = ConArea[k][1] = face_area[i];

        ConP1x[k] = face_normal_x[i];
        ConP1y[k] = face_normal_y[i];
        ConP1z[k] = face_normal_z[i];

        con_center_x[k] = polygon_center_x[i];
        con_center_y[k] = polygon_center_y[i];
        con_center_z[k] = polygon_center_z[i];
        con_normal_x[k] = face_normal_x[i];
        con_normal_y[k] = face_normal_y[i];
        con_normal_z[k] = face_normal_z[i];

        ConTr[k].resize(2);
        ConPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;
        NbTransmissibility++;
      }

      // if (ListV1[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M
      if (CodePolygon[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M
      {
        // std::cout << "M-F && F-M " << i << std::endl;

        connection_type[k] = ConnectionType::matrix_fracture;
        n_connection_elements[k] = 2;

        ConCV[k].resize(2);
        ConCV[k][0] = EQV[ListV1[i]];
        ConCV[k][1] = EQF[i];

        ConArea[k].resize(2);
        ConArea[k][0] = ConArea[k][1] = face_area[i];

        ConP1x[k] = face_normal_x[i];
        ConP1y[k] = face_normal_y[i];
        ConP1z[k] = face_normal_z[i];

        con_center_x[k] = polygon_center_x[i];
        con_center_y[k] = polygon_center_y[i];
        con_center_z[k] = polygon_center_z[i];
        con_normal_x[k] = face_normal_x[i];
        con_normal_y[k] = face_normal_y[i];
        con_normal_z[k] = face_normal_z[i];

        ConTr[k] .resize(2);
        ConPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;

        connection_type[k] = ConnectionType::matrix_fracture;
        n_connection_elements[k] = 2;
        ConCV[k].resize(2);
        ConCV[k][0] = EQV[ListV2[i]];
        ConCV[k][1] = EQF[i];

        ConArea[k].resize(2);
        ConArea[k][0] = ConArea[k][1] = face_area[i];

        ConP1x[k] = face_normal_x[i];
        ConP1y[k] = face_normal_y[i];
        ConP1z[k] = face_normal_z[i];

        con_center_x[k] = polygon_center_x[i];
        con_center_y[k] = polygon_center_y[i];
        con_center_z[k] = polygon_center_z[i];
        con_normal_x[k] = face_normal_x[i];
        con_normal_y[k] = face_normal_y[i];
        con_normal_z[k] = face_normal_z[i];

        ConTr[k] .resize(2);
        ConPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);
        k++;
        NbTransmissibility += 2;
      }
    }

    // std::cout << "F-F" << std::endl;
    // F-F Connections
    std::size_t i=0;
    while ( static_cast<int>(i) < NbEdges-1 )
    {
      j = i+1;
      while ( ( j<NbEdges ) && ( ListE1[i]==ListE1[j] ) && ( ListE2[i]==ListE2[j] ) )
        j++;

      if ( ( j-i ) >= 2 )  // otherwise it's the same face i think
      {
        connection_type[k] = ConnectionType::fracture_fracture;
        n_connection_elements[k] = ( j-i );

        ConCV[k].resize(j-i);
        for ( std::size_t n=i; n<j; n++ )  // CVs are faces here
          ConCV[k][n-i] = EQF[ListF[n]];

        con_center_x[k] = X[ListE1[i]];
        con_center_y[k] = Y[ListE1[i]];
        con_center_z[k] = Z[ListE1[i]];
        // i have no idea, it's not pependicular to the edge is it
        con_normal_x[k] = X[ListE1[i]] - X[ListE2[i]];
        con_normal_y[k] = Y[ListE1[i]] - Y[ListE2[i]];
        con_normal_z[k] = Z[ListE1[i]] - Z[ListE2[i]];

        ConArea[k].resize(j-1);

        for ( std::size_t n = i; n < j; n++ ) // Double check the formula
        {
          // what the damn hell???
          ConArea[k][n - i] = ZVolumeFactor[CVZone[ConCV[k][n - i]]] *
              sqrt ( con_normal_x[k] * con_normal_x[k] +
                     con_normal_y[k] * con_normal_y[k] +
                     con_normal_z[k] * con_normal_z[k] );
          // TODO TIMUR (F-F connection)
          ConArea[k][n - i] *= vTimurConnectionFactor[CVZone[ConCV[k][n - i]]];
        }
        // ConTr[k] = ( double* ) malloc ( ( j-i ) * ( sizeof ( double ) ) );
        ConTr[k].resize(j-i);
        ConPerm[k].resize(j-i);

        ConGeom[k].resize(j-i);
        ConMult[k].resize(j-i);

        NbTransmissibility += ( n_connection_elements[k] * ( n_connection_elements[k]-1 ) ) /2;
        k++;
      }
      i = j;
    }


}
/********************************************************************/
void CalcTranses::VolumeCorrection()  // Volume should be at least twice bigger...
{
    for (std::size_t i=0; i<NbPolygons; i++)
    {
      // if (CodePolygon[i] >= 0 && ListV2[i] < 0) // M-F ///
      if (ListV1[i] >= 0 && ListV2[i] < 0) // M-F ///
      {
        if (CVVolume[EQV[ListV1[i]]] > 2*CVVolume[EQF[i]])
          CVVolume[EQV[ListV1[i]]] -= CVVolume[EQF[i]];
      }

      // if (CodePolygon[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M ///
      if (ListV1[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M ///
      {
        if (CVVolume[EQV[ListV1[i]]] > 2*.5*CVVolume[EQF[i]])
          CVVolume[EQV[ListV1[i]]] -= .5*CVVolume[EQF[i]];

        if (CVVolume[EQV[ListV2[i]]] > 2*.5*CVVolume[EQF[i]])
          CVVolume[EQV[ListV2[i]]] -= .5*CVVolume[EQF[i]];
      }
    }

}
/********************************************************************/
void CalcTranses::ComputeContinuityNode()
{
    int j, k;
    double  hx,hy,hz,px,py,pz;

    con_center_proj_x.resize(NbConnections);
    con_center_proj_y.resize(NbConnections);
    con_center_proj_z.resize(NbConnections);

    for (std::size_t i=0; i<NbConnections; i++)
    {

      con_center_proj_x[i] = con_center_proj_y[i] = con_center_proj_z[i] = 0;

      if (connection_type[i] == ConnectionType::matrix_matrix)
      {
        px = CV_center_x[ConCV[i][1]] - CV_center_x[ConCV[i][0]];
        py = CV_center_y[ConCV[i][1]] - CV_center_y[ConCV[i][0]];
        pz = CV_center_z[ConCV[i][1]] - CV_center_z[ConCV[i][0]];

        // get projection of connection center onto line connecting two CVs
        ProjectionA(CV_center_x[ConCV[i][0]], CV_center_y[ConCV[i][0]], CV_center_z[ConCV[i][0]],
                    px, py, pz,
                    con_center_x[i], con_center_y[i], con_center_z[i],
                    con_normal_x[i], con_normal_y[i], con_normal_z[i],
                    &con_center_proj_x[i], &con_center_proj_y[i], &con_center_proj_z[i]);
      }
      else if (connection_type[i] == ConnectionType::matrix_fracture)
      {
        // get projection of connection center onto line connecting two CVs
        // in this case it's simply the center of the fracture
        con_center_proj_x[i] = CV_center_x[ConCV[i][1]];
        con_center_proj_y[i] = CV_center_y[ConCV[i][1]];
        con_center_proj_z[i] = CV_center_z[ConCV[i][1]];
      }
      else if (connection_type[i] == ConnectionType::fracture_fracture)
      {
        for (j=0; j<n_connection_elements[i]; j++)
        {
          // in this case this projection is somewhat compliated
          // probably it's a part of the star transformation?
          ProjectionB( CV_center_x[ConCV[i][j]], CV_center_y[ConCV[i][j]], CV_center_z[ConCV[i][j]],
                       con_center_x[i], con_center_y[i], con_center_z[i],
                       con_normal_x[i], con_normal_y[i], con_normal_z[i],
                       &hx, &hy, &hz);
          con_center_proj_x[i] += hx;
          con_center_proj_y[i] += hy;
          con_center_proj_z[i] += hz;
        }

        con_center_proj_x[i] = con_center_proj_x[i] / n_connection_elements[i];
        con_center_proj_y[i] = con_center_proj_y[i] / n_connection_elements[i];
        con_center_proj_z[i] = con_center_proj_z[i] / n_connection_elements[i];
      }
    }
}
/********************************************************************/
void CalcTranses::ComputeDirectionalPermeability()
{
    int j,k;
    double  fx,fy,fz,fl;

    for (std::size_t i=0; i<NbConnections; i++)
    {
      if (connection_type[i] == ConnectionType::matrix_matrix)
      {
        k = ConCV[i][0];
        assert(k >= 0);

        fx = CV_center_x[k] - con_center_proj_x[i];
        fy = CV_center_y[k] - con_center_proj_y[i];
        fz = CV_center_z[k] - con_center_proj_z[i];

        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        Kx = ZPermeability[CVZone[k]][0]*fx +
             ZPermeability[CVZone[k]][3]*fy +
             ZPermeability[CVZone[k]][4]*fz;

        Ky = ZPermeability[CVZone[k]][3]*fx+
             ZPermeability[CVZone[k]][1]*fy +
             ZPermeability[CVZone[k]][5]*fz;

        Kz = ZPermeability[CVZone[k]][4]*fx +
             ZPermeability[CVZone[k]][5]*fy +
             ZPermeability[CVZone[k]][2]*fz;

        ConPerm[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        // @HACK
        // We dont need a REAL value of conductivity
        // We just calculate geometric part (@HACK)
        Kx=1.0; Ky=1.0; Kz=1.0;
        if(ZConduction[CVZone[k]][0] != 0.0)
        {
          Kx = (ZConduction[CVZone[k]][0]*fx +
                ZConduction[CVZone[k]][3]*fy +
                ZConduction[CVZone[k]][4]*fz) / ZConduction[CVZone[k]][0];

          Ky = (ZConduction[CVZone[k]][3]*fx +
                ZConduction[CVZone[k]][1]*fy +
                ZConduction[CVZone[k]][5]*fz) / ZConduction[CVZone[k]][0];

          Kz = (ZConduction[CVZone[k]][4]*fx +
                ZConduction[CVZone[k]][5]*fy +
                ZConduction[CVZone[k]][2]*fz) / ZConduction[CVZone[k]][0];
        }

        ConMult[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        k = ConCV[i][1];
        assert(k >= 0);

        fx = CV_center_x[k] - con_center_proj_x[i];
        fy = CV_center_y[k] - con_center_proj_y[i];
        fz = CV_center_z[k] - con_center_proj_z[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        // WTF indices????
        // Ok it's not voigt notation
        // [0-2] are xyz
        // [3-5] are symmetric components
        Kx = ZPermeability[CVZone[k]][0] * fx +
             ZPermeability[CVZone[k]][3] * fy +
             ZPermeability[CVZone[k]][4] * fz;

        Ky = ZPermeability[CVZone[k]][3] * fx +
             ZPermeability[CVZone[k]][1] * fy +
             ZPermeability[CVZone[k]][5] * fz;

        Kz = ZPermeability[CVZone[k]][4] * fx +
             ZPermeability[CVZone[k]][5] * fy +
             ZPermeability[CVZone[k]][2] * fz;

        ConPerm[i][1] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        // std::cout << ConCV[i][0] << " " << ConCV[i][1] << "\t"
        //           << ConPerm[i][0] << " " << ConPerm[i][0] << std::endl;

        // @HACK
        // We dont need a REAL value of conductivity
        // We just calculate geometric part
        Kx = 1.0;
        Ky = 1.0;
        Kz = 1.0;

        if(ZConduction[CVZone[k]][0] != 0.0)
        {
          Kx = (ZConduction[CVZone[k]][0]*fx +
                ZConduction[CVZone[k]][3]*fy +
                ZConduction[CVZone[k]][4]*fz) / ZConduction[CVZone[k]][0];

          Ky = (ZConduction[CVZone[k]][3]*fx +
                ZConduction[CVZone[k]][1]*fy +
                ZConduction[CVZone[k]][5]*fz) / ZConduction[CVZone[k]][0];

          Kz = (ZConduction[CVZone[k]][4]*fx +
                ZConduction[CVZone[k]][5]*fy +
                ZConduction[CVZone[k]][2]*fz) / ZConduction[CVZone[k]][0];
        }
        ConMult[i][1]=sqrt(Kx*Kx+Ky*Ky+Kz*Kz);
      }
      if (connection_type[i] == ConnectionType::matrix_fracture)
      {
        // std::cout << "M-F" << std::endl;

        k = ConCV[i][0];
        assert(k >= 0);

        fx = CV_center_x[k] - con_center_proj_x[i];
        fy = CV_center_y[k] - con_center_proj_y[i];
        fz = CV_center_z[k] - con_center_proj_z[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        Kx = ZPermeability[CVZone[k]][0]*fx +
             ZPermeability[CVZone[k]][3]*fy +
             ZPermeability[CVZone[k]][4]*fz;

        Ky = ZPermeability[CVZone[k]][3]*fx +
             ZPermeability[CVZone[k]][1]*fy +
             ZPermeability[CVZone[k]][5]*fz;

        Kz = ZPermeability[CVZone[k]][4]*fx +
             ZPermeability[CVZone[k]][5]*fy +
             ZPermeability[CVZone[k]][2]*fz;

        ConPerm[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        k = ConCV[i][1];
        assert(k >= 0);

        ConPerm[i][1] = ZPermeability[CVZone[k]][0];    // Kn (0)

        // @HACK
        // We dont need a REAL value of conductivity
        // We just calculate geometric part
        Kx=1.0; Ky=1.0; Kz=1.0;
        if(ZConduction[CVZone[k]][0] != 0.0)
        {
          Kx = (ZConduction[CVZone[k]][0]*fx +
                ZConduction[CVZone[k]][3]*fy +
                ZConduction[CVZone[k]][4]*fz) / ZConduction[CVZone[k]][0];

          Ky = (ZConduction[CVZone[k]][3]*fx +
                ZConduction[CVZone[k]][1]*fy +
                ZConduction[CVZone[k]][5]*fz) / ZConduction[CVZone[k]][0];

          Kz = (ZConduction[CVZone[k]][4]*fx +
                ZConduction[CVZone[k]][5]*fy +
                ZConduction[CVZone[k]][2]*fz) / ZConduction[CVZone[k]][0];
        }

        ConMult[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
        ConMult[i][1] = 1.0;
      }
      if (connection_type[i] == ConnectionType::fracture_fracture)
      {
        // std::cout << "F-F" << std::endl;
        for (std::size_t j=0; j<n_connection_elements[i]; j++)
        {
          k = ConCV[i][j];
          assert(k >= 0);

          ConPerm[i][j] = ZPermeability[CVZone[k]][0];  // Kp (1)

          // @HACK
          // for fractures all directions are identiacal
          // and multiplicator is UNIT
          ConMult[i][j] = 1.0;
        }
      }
    }
}
/********************************************************************/
void CalcTranses::ComputeTransmissibilityPart()
{
    // int i,j,k;
    double  nx,ny,nz,nl,fx,fy,fz,fl;

    for (std::size_t i=0; i<NbConnections; i++)
    {
      if (connection_type[i] == ConnectionType::matrix_matrix)
      {
        std::size_t k  = ConCV[i][0];
        // nx = ConP1x[i];
        // ny = ConP1y[i];
        // nz = ConP1z[i];

        fx = CV_center_x[k] - con_center_proj_x[i];
        fy = CV_center_y[k] - con_center_proj_y[i];
        fz = CV_center_z[k] - con_center_proj_z[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);

        // why? it's changed after
        // fx = fx/fl;
        // fy = fy/fl;
        // fz = fz/fl;

        //    ConTr[i][0]=ConArea[i][0]*ConPerm[i][0]*ABS(nx*fx+ny*fy+nz*fz)/fl;
        ConTr[i][0]   = ConArea[i][0] * ConPerm[i][0] / fl;
        ConGeom[i][0] = ConArea[i][0] * ConMult[i][0] / fl;

        k  = ConCV[i][1];
        // nx = ConP2x[i];
        // ny = ConP2y[i];
        // nz = ConP2z[i];

        fx = CV_center_x[k] - con_center_proj_x[i];
        fy = CV_center_y[k] - con_center_proj_y[i];
        fz = CV_center_z[k] - con_center_proj_z[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        // why?
        // fx = fx/fl;
        // fy = fy/fl;
        // fz = fz/fl;

        //    ConTr[i][1]=ConArea[i][1]*ConPerm[i][1]*ABS(nx*fx+ny*fy+nz*fz)/fl;
        ConTr[i][1]   = ConArea[i][1] * ConPerm[i][1] * 1./fl;
        ConGeom[i][1] = ConArea[i][1] * ConMult[i][1] * 1./fl;
      }
      else if (connection_type[i] == ConnectionType::matrix_fracture)
      {
        std::size_t k = ConCV[i][0];
        nx = ConP1x[i];
        ny = ConP1y[i];
        nz = ConP1z[i];

        fx = CV_center_x[k] - con_center_proj_x[i];
        fy = CV_center_y[k] - con_center_proj_y[i];
        fz = CV_center_z[k] - con_center_proj_z[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        ConTr[i][0]   = ConArea[i][0] * ConPerm[i][0] * 1./fl;
        ConGeom[i][0] = ConArea[i][0] * ConMult[i][0] * 1./fl;

        k = ConCV[i][1];

        // is it like half a fracture?????
        // units are different from cell tran!!!!!
        ConTr[i][1] = ConArea[i][1] * ConPerm[i][1] *
            1./(.5*ZVolumeFactor[CVZone[k]]);

        ConGeom[i][1] = ConArea[i][1] * ConMult[i][1] *
            1./(.5*ZVolumeFactor[CVZone[k]]);
      }
      else if (connection_type[i] == ConnectionType::fracture_fracture)
      {
        for (std::size_t j=0; j<n_connection_elements[i]; j++)
        {
          std::size_t k  = ConCV[i][j];
          fx = CV_center_x[k] - con_center_proj_x[i];
          fy = CV_center_y[k] - con_center_proj_y[i];
          fz = CV_center_z[k] - con_center_proj_z[i];
          fl = sqrt(fx*fx + fy*fy + fz*fz);
          fx = fx/fl;
          fy = fy/fl;
          fz = fz/fl;

          ConTr[i][j] = ConArea[i][j] * ConPerm[i][j] * 1./fl;
          ConGeom[i][j] = ConArea[i][j] * ConMult[i][j] * 1./fl;
          }
      }
    }
}
/********************************************************************/
void CalcTranses::ComputeTransmissibilityList()
{
  // std::cout <<  " calculate transmissibilities"  << std::endl;

  iTr.resize(NbTransmissibility);
  jTr.resize(NbTransmissibility);
  Tij.resize(NbTransmissibility);
  TConductionIJ.resize(NbTransmissibility);

    std::size_t k=0;
    for ( std::size_t i=0; i<NbConnections; i++ )
    {
      if ( connection_type[i] == ConnectionType::matrix_matrix )
      {
        iTr[k] = ConCV[i][0];
        jTr[k] = ConCV[i][1];

        if ( ConTr[i][0] + ConTr[i][1] != 0.0 )
        {
          Tij[k] = ( ConTr[i][0]*ConTr[i][1] ) / ( ConTr[i][0] + ConTr[i][1] );
        }

        if ( ConGeom[i][0] + ConGeom[i][1] != 0.0 )
          TConductionIJ[k] = ( ConGeom[i][0] * ConGeom[i][1] ) /
                             ( ConGeom[i][0] + ConGeom[i][1] );

        if ( TConductionIJ[k] < 0.0 )
        {
          cout << "M-M or M-F : Conduction is negative: check values" << endl;
          exit ( 0 );
        }
        k++;
      }
      if ( connection_type[i] == ConnectionType::matrix_fracture )
      {
        iTr[k] = ConCV[i][0];
        jTr[k] = ConCV[i][1];

        if ( ConTr[i][0] + ConTr[i][1] != 0.0 )
        {
          Tij[k] = ( ConTr[i][0]*ConTr[i][1] ) / ( ConTr[i][0] + ConTr[i][1] );
        }

        if ( ConGeom[i][0] + ConGeom[i][1] != 0.0 )
          TConductionIJ[k] = ( ConGeom[i][0] * ConGeom[i][1] ) /
                             ( ConGeom[i][0] + ConGeom[i][1] );
        // if ( ConGeom[i][0] + ConGeom[i][1] != 0.0 )
        //   TConductionIJ[k] = std::min( ConGeom[i][0], ConGeom[i][1] );

        if ( TConductionIJ[k] < 0.0 )
        {
          cout << "M-M or M-F : Conduction is negative: check values" << endl;
          exit ( 0 );
        }
        k++;
      }
      else if ( connection_type[i] == ConnectionType::fracture_fracture )
      {
        double SumTr = 0;
        double SumTr2 = 0;
        for ( std::size_t j=0; j<n_connection_elements[i]; j++ )
          SumTr+=ConTr[i][j];

        for ( std::size_t j=0; j<n_connection_elements[i]; j++ )
          SumTr2+=ConGeom[i][j];

        for ( std::size_t j=0; j<n_connection_elements[i]-1; j++ )
          for ( std::size_t n=j+1; n<n_connection_elements[i]; n++ )
          {
            iTr[k] = ConCV[i][j];
            jTr[k] = ConCV[i][n];
            Tij[k] = ( ConTr[i][j]*ConTr[i][n] ) / SumTr;

            TConductionIJ[k] = ( ConGeom[i][j]*ConGeom[i][n] ) / SumTr2;
            if ( TConductionIJ[k] < 0.0 )
            {
              cout << "M-M or M-F : Conduction is negative: check values" << endl;
              exit ( 0 );
            }
            k++;
          }
      }
      else
      {
        // cout << "Wrong connection type " << k << ": "  << connection_type[i]  << endl;
        // exit ( 0 );
      }

    }

}

/********************************************************************/
void CalcTranses::compute_flow_data()
{
    double m1x, m1y, m1z,
        p1x, p1y, p1z,
        m2x, m2y, m2z,
        p2x, p2y, p2z,
        ix, iy, iz,
        inx, iny, inz;

    double mx, my, mz, ux, uy, uz, vx, vy, vz, nx, ny, nz, xi, yi, zi, nl, h, fx, fy, fz, fl;
    int   iswap;
    double    areatmp, volumetmp;
    double    hx,hy,hz;
    int   NbActivePolygon;
    int   NbFeatureCode;

    EQF.resize(NbPolygons);
    CodePolygon.resize(NbPolygons);
    FNodes.resize(NbPolygons);

    NbEdges = 0;
    NbCVs = 0;  // number of control volumes?

    for (std::size_t i=0; i<NbPolygons; i++)
    {
      FNodes[i].resize(vvFNodes[i].size());

      // polygon node indices
      for (std::size_t j=0; j<vvFNodes[i].size(); j++)
        FNodes[i][j] = vvFNodes[i][j];

      CodePolygon[i] = vCodePolygon[i];

      if (CodePolygon[i] >= 0)
      {
        EQF[i] = NbCVs++;
        NbEdges += vvFNodes[i].size();
      }
      else EQF[i] = -1;
    }

    NbActivePolygon = NbCVs;

    // NbVFaces.resize(NbPolyhedra);
    vvVFaces.resize(NbPolyhedra);
    EQV.resize(NbPolyhedra);
    CodePolyhedron.resize(NbPolyhedra);

    for (std::size_t i=0; i<NbPolyhedra; i++)
    {
      CodePolyhedron[i] = vCodePolyhedron[i];
      if (CodePolyhedron[i] >= 0)
      {
        EQV[i] = NbCVs++;
      }
      else EQV[i] = -1;
    }

    ZoneCode.resize(NbZones);
    ZVolumeFactor.resize(NbZones);
    ZPorosity.resize(NbZones);
    ZPermCode.resize(NbZones);

    ZPermeability.resize(NbZones, std::vector<double>(6, 0));
    ZConduction.resize(NbZones, std::vector<double>(6, 0));

    for (std::size_t i=0; i<NbZones; i++)
    {
      // std::cout << "i = " << i << std::endl;
      ZoneCode[i] = vZoneCode[i];
      std::size_t j = ZoneCode[i];
      ZVolumeFactor[j] = vZVolumeFactor[j];

      ZPorosity[j] = vZPorosity[j];
      ZPermCode[j] = vZPermCode[j];

      ZPermeability[j][0] = vZPermeability[j * 3 + 0];
      ZPermeability[j][1] = vZPermeability[j * 3 + 1];
      ZPermeability[j][2] = vZPermeability[j * 3 + 2];

      ZPermeability[j][3] = 0.0;
      ZPermeability[j][4] = 0.0;
      ZPermeability[j][5] = 0.0;

      ZConduction[j][0] = vZConduction[j * 3 + 0];
      ZConduction[j][1] = vZConduction[j * 3 + 1];
      ZConduction[j][2] = vZConduction[j * 3 + 2];

      ZConduction[j][3] = 0.0;
      ZConduction[j][4] = 0.0;
      ZConduction[j][5] = 0.0;
   }

    ComputeBasicGeometry();
    ComputeControlVolumeList();
    PrepareConnectionList();
    ConstructConnectionList();

    if (NbOptions == 1) VolumeCorrection();

    // printf("NbCVs = %d\n",NbCVs);
    // printf("NbTransmissibility = %d\n",NbTransmissibility);

    ComputeContinuityNode();
    ComputeDirectionalPermeability();
    ComputeTransmissibilityPart();
    ComputeTransmissibilityList();

    //////////////////////////////////
    ///// Computing Total Volume /////
    //////////////////////////////////

    TotalVolume = 0;
    for (std::size_t i=0; i<NbCVs; i++)
      TotalVolume += CVVolume[i];

    // printf("TotalVolume=%e\n", TotalVolume);

    /* OUTPUT MAPPING FOR GEOMECHANICS  */

    // ASSUMING FIRST THE FEATURE CODES ARE NUMBERED
    // printf("Prepare reservoir-geomechanical mapping...\n");
    NbFeatureCode = 0;
    for ( std::size_t i = 0; i < NbPolygons; i++ )
    {
      if ( CodePolygon[i] > NbFeatureCode )
      {
        NbFeatureCode++;
      }
    }

    // Features
    vector<int> vPriority; vPriority.resize(NbPolygons);
    int kk = 1;
    for (std::size_t i = 0; i < NbPolygons; i++)
      if(CodePolygon[i] > -1)
      {
        vPriority[i] = kk;
      kk++;
      }
}


void CalcTranses::writeOutputFiles(const std::string & output_path) const
{
  /* Creates files:
   * fl_vol.txt
   * fl_poro.txt
   * fl_depth.txt
   * fl_tran.txt
   * fl_tran_n.txt
   * fl_tran_u.txt
   */

  stringstream out;
  out << "model/";
  out << "/";
  string outputPath_ = output_path;// out.str();
 //////////////////////////
  ///// OUTPUT Volumes /////
  //////////////////////////
  int i;
  printf ( "Output Volumes...\n" );
  string outfile = outputPath_ + "fl_vol.txt";
  FILE * poutfile;
  poutfile=fopen ( outfile.c_str(),"w");
  //fprintf(out,"%d\n",NbCVs);
  fprintf ( poutfile,"%s\n","VOLUME" );
  for ( i=0; i<NbCVs; i++ )
    fprintf ( poutfile,"%e\n", CVVolume[i] );
  fprintf(poutfile,"%s\n","/\n");
  fclose(poutfile);
  ///////////////////////////
  ///// OUTPUT Porosity /////
  ///////////////////////////

    printf("Output Porosity...\n");
    outfile = outputPath_ + "fl_poro.txt";
    poutfile=fopen( outfile.c_str(),"w");
    fprintf(poutfile,"%s\n","PORO");
    for (i=0;i<NbCVs;i++) fprintf(poutfile,"%e\n",ZPorosity[CVZone[i]]);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

  ///////////////////////////
  ///// OUTPUT Depth  /////
  ///////////////////////////
    printf("Output Depth (z)...\n");
    outfile = outputPath_ + "fl_depth.txt";
    poutfile=fopen( outfile.c_str(),"w");
    //fprintf(out,"%d\n",NbCVs);
    fprintf(poutfile,"%s\n","DEPTH");
    for (i=0;i<NbCVs;i++)
      fprintf(poutfile,"%e\n",-CV_center_z[i]);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

    /* OUTPUT Transmissibility */
    printf("Output Transmissibility...\n");
    outfile = outputPath_ + "fl_tran.txt";
    poutfile=fopen( outfile.c_str(),"w");
    fprintf(poutfile,"%s\n","TPFACONNS");
    fprintf(poutfile,"%zu\n",NbTransmissibility);

    for (i=0;i<NbTransmissibility;i++)
      fprintf(poutfile,
              "%d\t%d\t%e\n",
              iTr[i],jTr[i],Tij[i] * transmissibility_conversion_factor);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

    printf("Output Transmissibility and Geometrical part...\n");
    outfile = outputPath_ + "fl_tran_n.txt";
    poutfile=fopen( outfile.c_str(),"w");
    fprintf(poutfile,"%s\n","TPFACONNSN");
    fprintf(poutfile,"%zu\n",NbTransmissibility);

    for (i=0;i<NbTransmissibility;i++)
      fprintf(poutfile,
              "%d\t%d\t%e\t%e\n",
              iTr[i],jTr[i],Tij[i] * transmissibility_conversion_factor,
              TConductionIJ[i]);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

    string outstring =  outputPath_ + "fl_tran_u.txt";
    poutfile = fopen(outstring.c_str(),"w");
    fprintf(poutfile,"GMUPDATETRANS\n");

    int k=0;
    for(i=0;i<NbConnections;i++)
    {
      if(connection_type[i]==1 || connection_type[i]==2)  // M-M, M-F /
      {
        if(connection_type[i]==1)   // M-M ////////////////////////////
        {
          // Con# connection_type i ai j aj -> Tij=ai*aj/(ai+aj)
          fprintf(poutfile,"%d\t%d\t%d\t%e\t%d\t%e\n", k,
                  connection_type[i], ConCV[i][0], ConTr[i][0], ConCV[i][1], ConTr[i][1]);
        }
        if(connection_type[i]==2)   // M-F ////////////////////////////
        {
          //Con# connection_type m am i ci ei ki ai=ci*ki/ei-> Tmi=am*ai/(am+ai)
        fprintf(poutfile,"%d\t%d\t%d\t%e\t%d\t%e\t%e\t%e\n", k,
                connection_type[i], ConCV[i][0], ConTr[i][0], ConCV[i][1],
                2.*ConArea[i][1], ZVolumeFactor[CVZone[ConCV[i][1]]], ConPerm[i][1]);
        }
      k++;
      }
      if(connection_type[i]==3) // F-F /////////////////////////////////////////////////
      {
        //Con# connection_type i ci ei ki j cj ej kj N n cn en kn
        //ai=ci*ki*ei aj=cj*kj*ej an=cn*kn*en
        //-> Tij=ai*aj/(SUM an)
        //
        for(int j=0;j<n_connection_elements[i]-1;j++)
          for(int n=j+1;n<n_connection_elements[i];n++)
          {
            fprintf(poutfile,"%d\t%d\t%d\t%e\t%e\t%e\t%d\t%e\t%e\t%e\t",
                    k,
                    connection_type[i],

                    ConCV[i][j],
                    ConTr[i][j]/(ConPerm[i][j]*ZVolumeFactor[CVZone[ConCV[i][j]]]),
                    ZVolumeFactor[CVZone[ConCV[i][j]]],
                    ConPerm[i][j],

                    ConCV[i][n],
                    ConTr[i][n]/(ConPerm[i][n]*ZVolumeFactor[CVZone[ConCV[i][n]]]),
                    ZVolumeFactor[CVZone[ConCV[i][n]]],
                    ConPerm[i][n]);

            fprintf(poutfile,"%d\t",n_connection_elements[i]);
            for(int m=0;m<n_connection_elements[i];m++)
              fprintf(poutfile,"%d\t%e\t%e\t%e\t",
                      ConCV[i][m],
                      ConTr[i][m]/(ConPerm[i][m]*ZVolumeFactor[CVZone[ConCV[i][m]]]),
                      ZVolumeFactor[CVZone[ConCV[i][m]]],
                      ConPerm[i][m]);

            fprintf(poutfile,"\n");
            k++;
          }
    }
    }
    fprintf(poutfile,"/\n");
    fclose(poutfile);
}


void CalcTranses::extractData(FlowData & data) const
{
  // std::cout << "extracting volume data" << std::endl;
  // Extract Volumes, porosity, depth
  data.cells.resize(NbCVs);
  // data.volumes.resize(NbCVs);
  // data.poro.resize(NbCVs);
  // data.depth.resize(NbCVs);
  for (std::size_t i=0; i<NbCVs; i++ )
  {
    data.cells[i].volume   = CVVolume[i];
    data.cells[i].porosity = ZPorosity[CVZone[i]];
    data.cells[i].depth    = -CV_center_z[i];
  }

  // Transmissibility
  for (std::size_t i=0;i<NbTransmissibility; i++)
  {
    auto & face = data.insert_connection(iTr[i], jTr[i]);
    face.transmissibility = Tij[i];
    face.thermal_conductivity = TConductionIJ[i];
  }

  // Geomechanics
  // data.connection_type.resize(NbConnections);
  // for(std::size_t i=0;i<NbConnections;i++)
  // {
  //   // data.connection_type[i] = connection_type[i];
  // }
}


void CalcTranses::save_output(const FlowData    & data,
                              const std::string & output_dir)
{
  const std::string fname_cell_data = "fl_cell_data.txt";
  const std::string fname_face_data = "fl_face_data.txt";

  {  // Write cell data
    ofstream out;
    out.open((output_dir + fname_cell_data).c_str());

    ///// OUTPUT Dimensions /////
    out << "DIMENS" << std::endl;
    out << data.cells.size() << "\t"
        << 1 << "\t" << 1 << "\t"
        << std::endl;
    out << "/" << std::endl << std::endl;

    ///// OUTPUT Volumes /////
    out << "VOLUME" << std::endl;
    for (const auto & v : data.cells)
      out << v.volume << std::endl;
    out << "/" << std::endl << std::endl;

    ///// OUTPUT Porosity /////
    out << "PORO" << std::endl;
    for (const auto & v : data.cells)
      out << v.porosity << std::endl;
    out << "/" << std::endl << std::endl;

    ///// OUTPUT Depth  /////
    out << "DEPTH" << std::endl;
    for (const auto & v : data.cells)
      out << v.depth << std::endl;
    out << "/" << std::endl << std::endl;

    // additional data (if any)
    for (std::size_t i=0; i<data.custom_names.size(); ++i)
    {
      out << data.custom_names[i] << std::endl;
      for (const auto & cell : data.cells)
      {
        if (cell.custom.size() != data.custom_names.size())
          assert(cell.custom.size() == data.custom_names.size());
        out << cell.custom[i] << std::endl;
      }
      out << "/" << std::endl << std::endl;
    }

    out.close();
  }

  { // Face data

    ofstream out;
    out.open((output_dir + fname_face_data).c_str());

    /* OUTPUT Transmissibility */
    out << "TPFACONNS" << std::endl;
    std::size_t n_connections = data.map_connection.size();
    out << n_connections << std::endl;
    for (const auto & conn : data.map_connection)
    {
      // const std::size_t iconn = conn.second;
      const auto element_pair = data.invert_hash(conn.first);
      const auto face = conn.second;
      out << element_pair.first << "\t"
          << element_pair.second << "\t"
          << std::scientific
          << face.transmissibility * transmissibility_conversion_factor
          << std::defaultfloat << std::endl;
    }
    out << "/" << std::endl;

    out.close();
  }

}

}
