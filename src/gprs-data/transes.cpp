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
  vzVolumeFactor.resize(NbZones);
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
void CalcTranses::ProjectionA( double mx,double my,double mz,
                               double px,double py,double pz,
                               double ix,double iy,double iz,
                               double nx,double ny,double nz,
                               double *hx,double *hy,double *hz)
{
  double  t;

  t = ((ix-mx)*nx + (iy-my)*ny + (iz-mz)*nz) / (px*nx + py*ny + pz*nz);
  *hx = mx + t*px;
  *hy = my + t*py;
  *hz = mz + t*pz;
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
    int k;
    double  xi, yi, zi, areatmp, volumetmp;
    double  ux, uy, uz,
        vx, vy, vz,
        nx,ny,nz,
        nl,h;

//////////////////////////////////////////////////////////////////
///// Polygon Area, Center of Mass, and Normal (Unit vector) /////
//////////////////////////////////////////////////////////////////

    FArea.resize(NbPolygons);
    FXG.resize(NbPolygons);
    FYG.resize(NbPolygons);
    FZG.resize(NbPolygons);
    Fnx.resize(NbPolygons);
    Fny.resize(NbPolygons);
    Fnz.resize(NbPolygons);

    for (std::size_t i=0; i<NbPolygons; i++)
    {
      FArea[i] = FXG[i] = FYG[i] = FZG[i] = Fnx[i] = Fny[i] = Fnz[i] = 0;
      for (std::size_t j=1; j<vvFNodes[i].size()-1; j++)
      {
        ux = X[FNodes[i][j]] - X[FNodes[i][0]];
        uy = Y[FNodes[i][j]] - Y[FNodes[i][0]];
        uz = Z[FNodes[i][j]] - Z[FNodes[i][0]];
        vx = X[FNodes[i][j+1]] - X[FNodes[i][0]];
        vy = Y[FNodes[i][j+1]] - Y[FNodes[i][0]];
        vz = Z[FNodes[i][j+1]] - Z[FNodes[i][0]];
        nx = (uy*vz - vy*uz);
        ny = (vx*uz - ux*vz);
        nz = (ux*vy - uy*vx);
        areatmp = .5*sqrt(nx*nx + ny*ny + nz*nz);

        FArea[i] += areatmp;
        FXG[i] += areatmp*(X[FNodes[i][0]] + X[FNodes[i][j]] + X[FNodes[i][j+1]])/3.;
        FYG[i] += areatmp*(Y[FNodes[i][0]] + Y[FNodes[i][j]] + Y[FNodes[i][j+1]])/3.;
        FZG[i] += areatmp*(Z[FNodes[i][0]] + Z[FNodes[i][j]] + Z[FNodes[i][j+1]])/3.;

        Fnx[i] +=.5*nx;
        Fny[i] +=.5*ny;
        Fnz[i] +=.5*nz;
      }

      FXG[i] = FXG[i] / FArea[i];
      FYG[i] = FYG[i] / FArea[i];
      FZG[i] = FZG[i] / FArea[i];

      Fnx[i] = Fnx[i] / FArea[i];
      Fny[i] = Fny[i] / FArea[i];
      Fnz[i] = Fnz[i] / FArea[i];

      nl = sqrt(Fnx[i]*Fnx[i] + Fny[i]*Fny[i] + Fnz[i]*Fnz[i]);

      Fnx[i] = Fnx[i] / nl;
      Fny[i] = Fny[i] / nl;
      Fnz[i] = Fnz[i] / nl;
    }

////////////////////////////////////////////////
///// Polyhedron Volume and Center of Mass /////
////////////////////////////////////////////////

    VVolume.resize(NbPolyhedra);
    VXG.resize(NbPolyhedra);
    VYG.resize(NbPolyhedra);
    VZG.resize(NbPolyhedra);

    for (std::size_t i=0; i<NbPolyhedra; i++)
    {
        VVolume[i] = VXG[i] = VYG[i] = VZG[i] = 0;
        xi = yi = zi = 0;

        for (std::size_t j=0; j<vvVFaces[i].size(); j++) // Defining a node inside the polyhedron
        {
            xi += FXG[vvVFaces[i][j]];
            yi += FYG[vvVFaces[i][j]];
            zi += FZG[vvVFaces[i][j]];
        }

        const double n_element_nodes = vvVFaces[i].size();
        xi = xi / n_element_nodes;
        yi = yi / n_element_nodes;
        zi = zi / n_element_nodes;

        for (std::size_t j=0; j<vvVFaces[i].size(); j++)
        {
          k = vvVFaces[i][j];
          h = Fnx[k] * (FXG[k]-xi) +
              Fny[k] * (FYG[k]-yi) +
              Fnz[k] * (FZG[k]-zi);

          volumetmp = fabs(h*FArea[k]) / 3.;

          if (std::isnan(volumetmp))
          {
            std::cout << "found nan volume" << std::endl;
            std::cout << "Polyhedron = " << i << std::endl;
            std::cout << "k = " << k << std::endl;
            std::cout << "h = " << h << std::endl;
            std::cout << "xi = " << xi << std::endl;
            std::cout << "yi = " << yi << std::endl;
            std::cout << "zi = " << zi << std::endl;
            std::cout << "n_element_nodes = " << n_element_nodes << std::endl;
            for (std::size_t j=0; j<vvVFaces[i].size(); j++) // Defining a node inside the polyhedron
            {
              std::cout << "FXG[vvVFaces[i][j]]; = " << FXG[vvVFaces[i][j]] << std::endl;
              std::cout << "FYG[vvVFaces[i][j]]; = " << FYG[vvVFaces[i][j]] << std::endl;
              std::cout << "FzG[vvVFaces[i][j]]; = " << FZG[vvVFaces[i][j]] << std::endl;
            }
            abort();
          }

          VXG[i] += (FXG[k] + .25*(xi-FXG[k])) * volumetmp;
          VYG[i] += (FYG[k] + .25*(yi-FYG[k])) * volumetmp;
          VZG[i] += (FZG[k] + .25*(zi-FZG[k])) * volumetmp;
          VVolume[i] += volumetmp;
        }

        VXG[i] = VXG[i] / VVolume[i];
        VYG[i] = VYG[i] / VVolume[i];
        VZG[i] = VZG[i] / VVolume[i];

    }
}
/********************************************************************/
void CalcTranses::ComputeControlVolumeList()
{
  // int j;

  CVType.resize(NbCVs);
  CVZone.resize(NbCVs);
  CVVolume.resize(NbCVs);
  CVx.resize(NbCVs);
  CVy.resize(NbCVs);
  CVz.resize(NbCVs);

  for (std::size_t i=0; i<NbPolygons; i++)
    if (EQF[i] != -1) // Active polygon
    {
      std::size_t j = EQF[i];
      CVZone[j] = CodePolygon[i];
      CVType[j] = 2;  // Feature
      CVx[j] = FXG[i];
      CVy[j] = FYG[i];
      CVz[j] = FZG[i];
      CVVolume[j] = FArea[i] * zVolumeFactor[CVZone[j]];
    }

  for (int i=0;i<NbPolyhedra;i++)
  {
    if (EQV[i] != -1) // Active polyhedron
    {
      std::size_t j = EQV[i];
      CVZone[j] = CodePolyhedron[i];
      CVType[j] = 1;  // Volume
      CVx[j] = VXG[i];
      CVy[j] = VYG[i];
      CVz[j] = VZG[i];
      CVVolume[j] = VVolume[i] * zVolumeFactor[CVZone[j]];
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
    if (NbPolyhedra > 0) // syncronize this part as the same as in the ConstructConnectionList
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

    conType.resize(NbConnections);
    conCV.resize(NbConnections);
    conTr.resize(NbConnections);
    conArea.resize(NbConnections);

    ConGeom.resize(NbConnections);
    ConMult.resize(NbConnections);

    conPerm.resize(NbConnections);
    ConP1x.resize(NbConnections);
    ConP1y.resize(NbConnections);
    ConP1z.resize(NbConnections);

    ConP2x.resize(NbConnections);
    ConP2y.resize(NbConnections);
    ConP2z.resize(NbConnections);

    ConIx.resize(NbConnections);
    ConIy.resize(NbConnections);
    ConIz.resize(NbConnections);

    ConVx.resize(NbConnections);
    ConVy.resize(NbConnections);
    ConVz.resize(NbConnections);


    NbTransmissibility = 0;

    k=0;
    // std::cout << "NbPolygons = " << NbPolygons << std::endl;

    if (NbPolyhedra > 0)
    for (std::size_t i=0; i<NbPolygons; i++)
    {
      if (CodePolygon[i] < 0 && ListV2[i] >= 0)  // M-M
      {
        // std::cout << "M-M " << i << std::endl;
        conType[k] = 1;

        conCV[k].resize(2);
        conCV[k][0] = EQV[ListV1[i]];
        conCV[k][1] = EQV[ListV2[i]];

        conArea[k].resize(2);
        conArea[k][0] = conArea[k][1] = FArea[i];

        ConP1x[k] = Fnx[i];
        ConP1y[k] = Fny[i];
        ConP1z[k] = Fnz[i];

        ConP2x[k] = Fnx[i];
        ConP2y[k] = Fny[i];
        ConP2z[k] = Fnz[i];

        ConIx[k] = FXG[i];
        ConIy[k] = FYG[i];
        ConIz[k] = FZG[i];

        ConVx[k] = Fnx[i];
        ConVy[k] = Fny[i];
        ConVz[k] = Fnz[i];

        conTr[k].resize(2);
        conPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;
        NbTransmissibility++;
      }

      // if (ListV1[i] >= 0 && ListV2[i] < 0) // M-F
      if (CodePolygon[i] >= 0 && ListV2[i] < 0) // M-F
      {
        // std::cout << "M-F " << i << std::endl;
        conType[k] = 2;

        conCV[k].resize(2);

        if (ListV1[i] < 0)
        {
          std::cout << "polygon " << i << "\t";
          std::cout << "wrong ListV2 index " << ListV1[i] << std::endl;
          std::cout << "aborting" << std::endl;
          exit(0);
        }

        conCV[k][0] = EQV[ListV1[i]];
        conCV[k][1] = EQF[i];

        conArea[k].resize(2);
        conArea[k][0] = conArea[k][1] = FArea[i];

        ConP1x[k] = Fnx[i];
        ConP1y[k] = Fny[i];
        ConP1z[k] = Fnz[i];

        ConIx[k] = FXG[i];
        ConIy[k] = FYG[i];
        ConIz[k] = FZG[i];
        ConVx[k] = Fnx[i];
        ConVy[k] = Fny[i];
        ConVz[k] = Fnz[i];

        conTr[k].resize(2);
        conPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;
        NbTransmissibility++;
      }

      // if (ListV1[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M
      if (CodePolygon[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M
      {
        // std::cout << "M-F && F-M " << i << std::endl;

        conType[k] = 2;

        conCV[k].resize(2);
        conCV[k][0] = EQV[ListV1[i]];
        conCV[k][1] = EQF[i];

        conArea[k].resize(2);
        conArea[k][0] = conArea[k][1] = FArea[i];

        ConP1x[k] = Fnx[i];
        ConP1y[k] = Fny[i];
        ConP1z[k] = Fnz[i];

        ConIx[k] = FXG[i];
        ConIy[k] = FYG[i];
        ConIz[k] = FZG[i];
        ConVx[k] = Fnx[i];
        ConVy[k] = Fny[i];
        ConVz[k] = Fnz[i];

        conTr[k] .resize(2);
        conPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;

        conType[k] = 2;
        conCV[k].resize(2);
        conCV[k][0] = EQV[ListV2[i]];
        conCV[k][1] = EQF[i];

        conArea[k].resize(2);
        conArea[k][0] = conArea[k][1] = FArea[i];

        ConP1x[k] = Fnx[i];
        ConP1y[k] = Fny[i];
        ConP1z[k] = Fnz[i];

        ConIx[k] = FXG[i];
        ConIy[k] = FYG[i];
        ConIz[k] = FZG[i];
        ConVx[k] = Fnx[i];
        ConVy[k] = Fny[i];
        ConVz[k] = Fnz[i];

        conTr[k] .resize(2);
        conPerm[k].resize(2);

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
      {
        j++;
      }
      if ( ( j-i ) >= 2 )
      {
        conType[k] = 3;
        std::size_t conN = ( j-i );

        conCV[k].resize(j-i);
        for ( std::size_t n=i; n<j; n++ )
          conCV[k][n-i] = EQF[ListF[n]];

        ConIx[k] = X[ListE1[i]];
        ConIy[k] = Y[ListE1[i]];
        ConIz[k] = Z[ListE1[i]];
        ConVx[k] = X[ListE1[i]] - X[ListE2[i]];
        ConVy[k] = Y[ListE1[i]] - Y[ListE2[i]];
        ConVz[k] = Z[ListE1[i]] - Z[ListE2[i]];

        conArea[k].resize(j-1);

        for ( std::size_t n = i; n < j; n++ ) // Double check the formula
        {
          conArea[k][n - i] = zVolumeFactor[CVZone[conCV[k][n - i]]] *
              sqrt ( ConVx[k] * ConVx[k] + ConVy[k] * ConVy[k] + ConVz[k] * ConVz[k] );
          // TODO TIMUR (F-F connection)
          conArea[k][n - i] *= vTimurConnectionFactor[CVZone[conCV[k][n - i]]];
        }
        // conTr[k] = ( double* ) malloc ( ( j-i ) * ( sizeof ( double ) ) );
        conTr[k].resize(j-i);
        conPerm[k].resize(j-i);

        ConGeom[k].resize(j-i);
        ConMult[k].resize(j-i);

        NbTransmissibility += ( conN* ( conN-1 ) ) /2;
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

    Conhx.resize(NbConnections);
    Conhy.resize(NbConnections);
    Conhz.resize(NbConnections);

    for (std::size_t i=0; i<NbConnections; i++)
    {

      Conhx[i] = Conhy[i] = Conhz[i] = 0;

      if (conType[i] == 1)  // M-M ///
      {
        px = CVx[conCV[i][1]] - CVx[conCV[i][0]];
        py = CVy[conCV[i][1]] - CVy[conCV[i][0]];
        pz = CVz[conCV[i][1]] - CVz[conCV[i][0]];

        ProjectionA(  CVx[conCV[i][0]], CVy[conCV[i][0]], CVz[conCV[i][0]],
                      px, py, pz,
                      ConIx[i], ConIy[i], ConIz[i],
                      ConVx[i], ConVy[i], ConVz[i],
                      &Conhx[i], &Conhy[i], &Conhz[i]);
      }
      else if (conType[i] == 2)  // M-F ///
      {
        Conhx[i] = CVx[conCV[i][1]];
        Conhy[i] = CVy[conCV[i][1]];
        Conhz[i] = CVz[conCV[i][1]];
      }
      else if (conType[i] == 3)  // F-F ///
      {
        std::size_t conN = conCV[i].size();
        for (j=0;j<conN;j++)
        {
          ProjectionB( CVx[conCV[i][j]], CVy[conCV[i][j]], CVz[conCV[i][j]],
                       ConIx[i], ConIy[i], ConIz[i],
                       ConVx[i], ConVy[i], ConVz[i],
                       &hx, &hy, &hz);
          Conhx[i] += hx;
          Conhy[i] += hy;
          Conhz[i] += hz;
        }

        Conhx[i] = Conhx[i] / conN;
        Conhy[i] = Conhy[i] / conN;
        Conhz[i] = Conhz[i] / conN;
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
      if (conType[i] == 1)  // M-M //
      {
        k = conCV[i][0];
        assert(k >= 0);

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        Kx = ZPermeability[CVZone[k]][0]*fx +
             ZPermeability[CVZone[k]][3]*fy+
             ZPermeability[CVZone[k]][4]*fz;

        Ky = ZPermeability[CVZone[k]][3]*fx+
             ZPermeability[CVZone[k]][1]*fy +
             ZPermeability[CVZone[k]][5]*fz;

        Kz = ZPermeability[CVZone[k]][4]*fx +
             ZPermeability[CVZone[k]][5]*fy +
             ZPermeability[CVZone[k]][2]*fz;

        conPerm[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

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

        k = conCV[i][1];
        assert(k >= 0);

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        Kx = ZPermeability[CVZone[k]][0] * fx +
             ZPermeability[CVZone[k]][3] * fy +
             ZPermeability[CVZone[k]][4] * fz;

        Ky = ZPermeability[CVZone[k]][3] * fx +
             ZPermeability[CVZone[k]][1] * fy +
             ZPermeability[CVZone[k]][5] * fz;

        Kz = ZPermeability[CVZone[k]][4] * fx +
             ZPermeability[CVZone[k]][5] * fy +
             ZPermeability[CVZone[k]][2] * fz;

        conPerm[i][1] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

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
      if (conType[i] == 2)  // M-F ////
      {
        // std::cout << "M-F" << std::endl;

        k = conCV[i][0];
        assert(k >= 0);

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
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

        conPerm[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        k = conCV[i][1];
        assert(k >= 0);

        conPerm[i][1] = ZPermeability[CVZone[k]][0];    // Kn (0)

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
      if (conType[i] == 3)  // F-F /////////////////////////////////////////////////
      {
        // std::cout << "F-F" << std::endl;
        std::size_t conN = conCV[i].size();
        for (std::size_t j=0; j<conN; j++)
        {
          k = conCV[i][j];
          assert(k >= 0);

          conPerm[i][j] = ZPermeability[CVZone[k]][0];  // Kp (1)

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
      if (conType[i] == 1)  // M-M ///
      {
        std::size_t k  = conCV[i][0];
        nx = ConP1x[i];
        ny = ConP1y[i];
        nz = ConP1z[i];

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        //    conTr[i][0]=conArea[i][0]*conPerm[i][0]*fabs(nx*fx+ny*fy+nz*fz)/fl;
        conTr[i][0]   = conArea[i][0]*conPerm[i][0] * 1./fl;
        ConGeom[i][0] = conArea[i][0]*ConMult[i][0] * 1./fl;

        k  = conCV[i][1];
        nx = ConP2x[i];
        ny = ConP2y[i];
        nz = ConP2z[i];

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        //    conTr[i][1]=conArea[i][1]*conPerm[i][1]*fabs(nx*fx+ny*fy+nz*fz)/fl;
        conTr[i][1]   = conArea[i][1] * conPerm[i][1] * 1./fl;
        ConGeom[i][1] = conArea[i][1] * ConMult[i][1] * 1./fl;
      }
      else if (conType[i] == 2)  // M-F //
      {
        std::size_t k = conCV[i][0];
        nx = ConP1x[i];
        ny = ConP1y[i];
        nz = ConP1z[i];

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        conTr[i][0]   = conArea[i][0] * conPerm[i][0] * 1./fl;
        ConGeom[i][0] = conArea[i][0] * ConMult[i][0] * 1./fl;

        k = conCV[i][1];

        conTr[i][1] = conArea[i][1] * conPerm[i][1] *
            1./(.5*zVolumeFactor[CVZone[k]]);

        ConGeom[i][1] = conArea[i][1] * ConMult[i][1] *
            1./(.5*zVolumeFactor[CVZone[k]]);
      }
      else if (conType[i] == 3)  // F-F //
      {
        std::size_t conN = conCV[i].size();
        for (std::size_t j=0; j<conN; j++)
        {
          std::size_t k  = conCV[i][j];
          fx = CVx[k] - Conhx[i];
          fy = CVy[k] - Conhy[i];
          fz = CVz[k] - Conhz[i];
          fl = sqrt(fx*fx + fy*fy + fz*fz);
          fx = fx/fl;
          fy = fy/fl;
          fz = fz/fl;

          conTr[i][j] = conArea[i][j] * conPerm[i][j] * 1./fl;
          ConGeom[i][j] = conArea[i][j] * ConMult[i][j] * 1./fl;
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
      if ( conType[i] == 1 ) // M-M //
      {
        iTr[k] = conCV[i][0];
        jTr[k] = conCV[i][1];

        if ( conTr[i][0] + conTr[i][1] != 0.0 )
        {
          Tij[k] = ( conTr[i][0]*conTr[i][1] ) / ( conTr[i][0] + conTr[i][1] );
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
      if ( conType[i] == 2 ) // M-F //
      {
        iTr[k] = conCV[i][0];
        jTr[k] = conCV[i][1];

        if ( conTr[i][0] + conTr[i][1] != 0.0 )
        {
          Tij[k] = ( conTr[i][0]*conTr[i][1] ) / ( conTr[i][0] + conTr[i][1] );
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
      else if ( conType[i] == 3 ) // F-F ///
      {
        double SumTr = 0;
        double SumTr2 = 0;
        std::size_t conN = conCV[i].size();
        for ( std::size_t j=0; j<conN; j++ )
          SumTr+=conTr[i][j];

        for ( std::size_t j=0; j<conN; j++ )
          SumTr2+=ConGeom[i][j];

        for ( std::size_t j=0; j<conN-1; j++ )
          for ( std::size_t n=j+1; n<conN; n++ )
          {
            iTr[k] = conCV[i][j];
            jTr[k] = conCV[i][n];
            Tij[k] = ( conTr[i][j]*conTr[i][n] ) / SumTr;

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
        // cout << "Wrong connection type " << k << ": "  << conType[i]  << endl;
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
    zVolumeFactor.resize(NbZones);
    ZPorosity.resize(NbZones);
    ZPermCode.resize(NbZones);

    ZPermeability.resize(NbZones, std::vector<double>(6, 0));
    ZConduction.resize(NbZones, std::vector<double>(6, 0));

    for (std::size_t i=0; i<NbZones; i++)
    {
      ZoneCode[i] = vZoneCode[i];
      std::size_t j = ZoneCode[i];
      zVolumeFactor[j] = vzVolumeFactor[j];

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
    fprintf(poutfile,"%s\n","DEPTH");
    for (i=0;i<NbCVs;i++)
      fprintf(poutfile,"%e\n",-CVz[i]);

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
      if(conType[i]==1 || conType[i]==2)  // M-M, M-F /
      {
        if(conType[i]==1)   // M-M ////////////////////////////
        {
          // Con# conType i ai j aj -> Tij=ai*aj/(ai+aj)
          fprintf(poutfile,"%d\t%d\t%d\t%e\t%d\t%e\n", k,
                  conType[i], conCV[i][0], conTr[i][0], conCV[i][1], conTr[i][1]);
        }
        if(conType[i]==2)   // M-F ////////////////////////////
        {
          //Con# conType m am i ci ei ki ai=ci*ki/ei-> Tmi=am*ai/(am+ai)
        fprintf(poutfile,"%d\t%d\t%d\t%e\t%d\t%e\t%e\t%e\n", k,
                conType[i], conCV[i][0], conTr[i][0], conCV[i][1],
                2.*conArea[i][1], zVolumeFactor[CVZone[conCV[i][1]]], conPerm[i][1]);
        }
      k++;
      }
      if(conType[i]==3) // F-F /////////////////////////////////////////////////
      {
        //Con# conType i ci ei ki j cj ej kj N n cn en kn
        //ai=ci*ki*ei aj=cj*kj*ej an=cn*kn*en
        //-> Tij=ai*aj/(SUM an)
        //
        std::size_t conN = conCV[i].size();
        for(int j=0;j<conN-1;j++)
          for(int n=j+1;n<conN;n++)
          {
            fprintf(poutfile,"%d\t%d\t%d\t%e\t%e\t%e\t%d\t%e\t%e\t%e\t",
                    k,
                    conType[i],

                    conCV[i][j],
                    conTr[i][j]/(conPerm[i][j]*zVolumeFactor[CVZone[conCV[i][j]]]),
                    zVolumeFactor[CVZone[conCV[i][j]]],
                    conPerm[i][j],

                    conCV[i][n],
                    conTr[i][n]/(conPerm[i][n]*zVolumeFactor[CVZone[conCV[i][n]]]),
                    zVolumeFactor[CVZone[conCV[i][n]]],
                    conPerm[i][n]);

            fprintf(poutfile,"%d\t",conN);
            for(int m=0;m<conN;m++)
              fprintf(poutfile,"%d\t%e\t%e\t%e\t",
                      conCV[i][m],
                      conTr[i][m]/(conPerm[i][m]*zVolumeFactor[CVZone[conCV[i][m]]]),
                      zVolumeFactor[CVZone[conCV[i][m]]],
                      conPerm[i][m]);

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
  for (std::size_t i=0; i<NbCVs; i++ )
  {
    data.cells[i].volume   = CVVolume[i];
    data.cells[i].porosity = ZPorosity[CVZone[i]];
    data.cells[i].depth    = -CVz[i];
  }

  // Transmissibility
  for (std::size_t i=0;i<NbConnections; i++)
  {
    if(conType[i]==1 || conType[i]==2)  // M-M, M-F /
    {
        auto & face = data.insert_connection(iTr[i], jTr[i]);
        std::size_t conN = conCV[i].size();
        face.transmissibility = Tij[i];
        face.thermal_conductivity = TConductionIJ[i];
        face.conType = conType[i];
        face.conCV.resize(conN);
        face.conTr.resize(conN);
        face.conArea.resize(conN);
        face.conPerm.resize(conN);
        face.zVolumeFactor.resize(conN);
        for(std::size_t m=0; m < conN; m++){
            face.conCV[m] = conCV[i][m];
            face.conTr[m] = conTr[i][m];
            face.conArea[m] = conArea[i][m];
            face.conPerm[m] = conPerm[i][m];
            face.zVolumeFactor[m] = zVolumeFactor[CVZone[conCV[i][m]]];
        }
    } else{
        double SumTr = 0;
        double SumTr2 = 0;
        std::size_t conN = conCV[i].size();
        for ( std::size_t j=0; j<conN; j++ )
          SumTr+=conTr[i][j];

        for ( std::size_t j=0; j<conN; j++ )
          SumTr2+=ConGeom[i][j];

        for ( std::size_t j=0; j<conN-1; j++ ){
            for ( std::size_t n=j+1; n<conN; n++ )
            {// eventhough the name is face; treat as connection.
                // IN edfm intersection, acturally there in only one connection in frac. THINK ABOUT HOW TO APPROXIMATE THIS.
                std::size_t iTr_ = conCV[i][j];
                std::size_t jTr_ = conCV[i][n];
                double Tij_ = ( conTr[i][j]*conTr[i][n] ) / SumTr;
                double TConductionIJ_ = ( ConGeom[i][j]*ConGeom[i][n] ) / SumTr2;
                auto & face = data.insert_connection(iTr_, jTr_);
                face.transmissibility = Tij_;
                face.thermal_conductivity = TConductionIJ_;
                face.conType = conType[i];
                face.conCV.resize(conN);
                face.conTr.resize(conN);
                face.conArea.resize(conN);
                face.conPerm.resize(conN);
                face.zVolumeFactor.resize(conN);

                std::vector<std::size_t> vecIndex(conN);
                vecIndex[0] = j;
                vecIndex[1] = n;
                std::size_t tempI = 2;
                for(std::size_t m=0; m < conN; m++){
                    if(m!=j && m!=n){
                        vecIndex[tempI] = m;
                        tempI++;
                    }
                }

                for(std::size_t m=0; m < conN; m++){
                    face.conCV[m] = conCV[i][vecIndex[m]];
                    face.conTr[m] = conTr[i][vecIndex[m]];
                    face.conArea[m] = conArea[i][vecIndex[m]];
                    face.conPerm[m] = conPerm[i][vecIndex[m]];
                    face.zVolumeFactor[m] = zVolumeFactor[CVZone[conCV[i][vecIndex[m]]]];
                }
            }
        }
    }
  }
}


void CalcTranses::save_output(const FlowData    & data,
                              const std::string & output_dir)
{
  const std::string fname_cell_data = "fl_cell_data.txt";
  const std::string fname_face_data = "fl_face_data.txt";
  const std::string fname_gmupdate_data = "fname_trans_data.txt";

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

  { // Geometric part of transmissibilities
      ofstream out;
      out.open((output_dir + fname_gmupdate_data).c_str());


      out << "GMUPDATETRANS" << std::endl;
      int k=0;
      for (const auto & conn : data.map_connection)
        {
        const auto face = conn.second;
        if(face.conType==1 || face.conType==2)	// M-M, M-F ////////////////////////////
          {
          if(face.conType==1)		// M-M ////////////////////////////
          {
              // Con# conType i ai j aj -> Tij=ai*aj/(ai+aj)
              out << k << "\t"
                  << face.conType << "\t"
                  << face.conCV[0] << "\t"
                  << face.conTr[0] << "\t"
                  << face.conCV[1] << "\t"
                  << face.conTr[1] << std::endl;
          }
          if(face.conType==2)		// M-F ////////////////////////////
          {
              //Con# conType m am i ci ei ki ai=ci*ki/ei-> Tmi=am*ai/(am+ai)
              out << k << "\t"
                  << face.conType << "\t"
                  << face.conCV[0] << "\t"
                  << face.conTr[0] << "\t"
                  << face.conCV[1] << "\t"
                  << 2.*face.conArea[1]  << "\t"
                  << face.zVolumeFactor[1] << "\t"
                  << face.conPerm[1] << "\t"
                  << std::endl;
          }
          k++;
          }
        if(face.conType==3)	// F-F /////////////////////////////////////////////////
          {
           //Con# conType i ci ei ki j cj ej kj N n cn en kn
           //ai=ci*ki*ei aj=cj*kj*ej an=cn*kn*en
           //-> Tij=ai*aj/(SUM an)
           //
            std::size_t j = 0;
            std::size_t n = 1;
            std::size_t conN = face.conCV.size();
            out << k << "\t"
              << face.conType << "\t"
              << face.conCV[j] << "\t"
              << face.conTr[j]/(face.conPerm[j]*face.zVolumeFactor[j]) << "\t"
              << face.zVolumeFactor[j] << "\t"
              << face.conPerm[j]  << "\t"
              << face.conCV[n] << "\t"
              << face.conTr[n]/(face.conPerm[n]*face.zVolumeFactor[n])<< "\t"
              << face.zVolumeFactor[n] << "\t"
              << face.conPerm[n] << "\t"
              << face.conCV[n] << "\t";
            for(std::size_t m=0; m< conN; m++)
              out << face.conCV[m] << "\t"
                  << face.conTr[m]/(face.conPerm[m]*face.zVolumeFactor[m]) << "\t"
                  << face.zVolumeFactor[m] << "\t"
                  << face.conPerm[m] << "\t";
            out << std::endl;
            k++;
           }
        }
      out << "/" << std::endl;

      out.close();
  }

}

}
