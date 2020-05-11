
#include <algorithm>
#include <chrono>
#include "StabFEM.h"
#include "KimMoinFlow.h"
#include "UtilitiesGeneral.h"
//#include "SolutionData.h"
//#include "LagrangeElem2DNavierStokesTria3Node.h"
#include "LagrangeElem2DNavierStokesQuad4Node.h"
//#include "LagrangeElem3DNavierStokesTetra4Node.h"
//#include "LagrangeElem3DNavierStokesHexa8Node.h"
#include "FunctionsSolver.h"


using namespace std;


StabFEM::StabFEM()
{
    ndof = 0; nElem = 0; nNode = 0; npElem = 0; fileCount = 0;
    totalDOF = 0;

    AlgoType = 2;

    elems = nullptr;

}


StabFEM::~StabFEM()
{
    if(elems != NULL)
    {
      for(int ii=0;ii<nElem;++ii)
        delete elems[ii];

      delete [] elems;
      elems = NULL;
    }

  phase = -1; error = 0;
/*
  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &totalDOF, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);
*/
}




void  StabFEM::readInputData(string&  fname)
{
    cout << " Reading input data \n\n " << endl;

    infilename = fname;

    std::ifstream  infile( string(infilename+".dat") );

    if(infile.fail())
    {
        cout << " Could not open the input nodes file " << endl;
        exit(1);
    }

    string  line, stringVal, stringVec[10];
    int  ii, arrayInt[100];
    double  tempDbl;

    // read the dimension
    infile >> stringVal >> ndim;
    ndof = ndim+1;

    // read number of nodes per element
    infile >> stringVal >> npElem;

    // read number of nodes
    infile >> stringVal >> nNode;

    // read number of elements
    infile >> stringVal >> nElem;

    // read number of Dirichlet BCs
    infile >> stringVal >> nDBC;

    // read number of Force BCs
    infile >> stringVal >> nFBC;

    // read number of Output nodes
    infile >> stringVal >> nOutputFaceLoads;

    cout << " ndim              =  " << ndim << endl;
    cout << " nNode             =  " << nNode << endl;
    cout << " nElem             =  " << nElem << endl;
    cout << " npElem            =  " << npElem << endl;
    cout << " nDBC              =  " << nDBC << endl;
    cout << " nFBC              =  " << nFBC << endl;
    cout << " nOutputFaceLoads  =  " << nOutputFaceLoads << endl;

    // read nodal coordinates
    ////////////////////////////////////////////

    cout << " reading nodes " << endl;

    node_coords.resize(nNode);
    for(ii=0; ii<nNode; ++ii)
      node_coords[ii].resize(ndim);

    infile >> stringVal ;
    cout << " reading " << stringVal << endl;
    if(ndim == 2)
    {
      for(ii=0; ii<nNode; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] ;
      }
    }
    else
    {
      for(ii=0; ii<nNode; ++ii)
      {
        infile >> tempDbl >> node_coords[ii][0] >> node_coords[ii][1] >> node_coords[ii][2];
      }
    }

    // read elements
    ////////////////////////////////////////////
    infile >> stringVal ;
    cout << " reading " << stringVal << '\t' << npElem << endl;

    elemConn.resize(nElem);

    if(npElem == 3)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6];

        //printf("%6d \t %6d \t %6d \t %6d \n", arrayInt[4], arrayInt[5], arrayInt[6], arrayInt[7]);

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 4)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7];

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else if(npElem == 8)
    {
      for(int ee=0; ee<nElem; ++ee)
      {
        elemConn[ee].resize(npElem);

        infile >> arrayInt[0] >> arrayInt[1] >> arrayInt[2] >> arrayInt[3] >> arrayInt[4] >> arrayInt[5] >> arrayInt[6] >> arrayInt[7] >> arrayInt[8] >> arrayInt[9] >> arrayInt[10] >> arrayInt[11] ;

        for(ii=0; ii<npElem; ++ii)
          elemConn[ee][ii] = arrayInt[4+ii]-1;
      }
    }
    else
    {
      cerr << " Invalid npElem " << npElem << endl;
      exit(-1);
    }

    //
    // Read Dirichlet BC data
    //
    ////////////////////////////////////////////

    vector<double>  vecDblTemp(3);

    infile >> stringVec[0] >> stringVec[1] >> stringVec[2] ;
    cout << " reading " << stringVec[0] << stringVec[1] << stringVec[2] << endl;

    for(ii=0; ii<nDBC; ++ii)
    {
      infile >> arrayInt[0] >> arrayInt[1] >> tempDbl ;

      vecDblTemp[0] = arrayInt[0]-1;
      vecDblTemp[1] = arrayInt[1]-1;
      vecDblTemp[2] = tempDbl;

      DirichletBCs.push_back(vecDblTemp);
    }

    //
    // Read Output data
    //
    ////////////////////////////////////////////

    infile >> stringVal ;
    cout << " reading " << stringVal << endl;

    outputEdges.resize(nOutputFaceLoads);

    for(ii=0; ii<nOutputFaceLoads; ++ii)
    {
      outputEdges[ii].resize(1);

      //infile >> arrayInt[0] >> arrayInt[1];
      infile >> arrayInt[0];

      outputEdges[ii][0] = arrayInt[0]-1;
      //outputEdges[ii][1] = arrayInt[1]-1;
    }

    infile.close();


    fout_convdata.open(string(infilename+"-forces.dat"), ios::out | ios::trunc );

    if(fout_convdata.fail())
    {
       cout << " Could not open the Output file" << endl;
       exit(1);
    }

    fout_convdata.setf(ios::fixed);
    fout_convdata.setf(ios::showpoint);
    fout_convdata.precision(14);


    cout << " Input files have been read successfully \n\n " << endl;

    return;
}



void StabFEM::readControlParameters(string& fname)
{
    cout << " StabFEM::readControlParameters " << endl;

    //ifstream  infile("control-parameters.dat");
    ifstream  infile(fname);

    if(infile.fail())
    {
       cout << " Could not open 'control-parameters.dat' file " << endl;
       exit(-1);
    }

    //time integration parameters
    timeData[1] = 1.0;   timeData[2] = 0.0;

    string  stringVal;

    //density
    infile >> stringVal >> elemData[0];

    //viscosity
    infile >> stringVal >> elemData[1];

    //Body force in X-, Y- and Z- direction
    infile >> stringVal >> elemData[2] >> elemData[3] >> elemData[4];

    //Stabilisation: SUPG, PSPG, LSIC
    infile >> stringVal >> elemData[8] >> elemData[9] >> elemData[10];

    //tis
    infile >> stringVal >> tis;

    //rhoInf = 0.0;
    infile >> stringVal >> rhoInf;

    // timestep
    infile >> stringVal >> dt;

    // final time
    infile >> stringVal >> timeFinal;

    // Maximum number of time steps
    infile >> stringVal >> stepsMax;

    // Output file frequency
    infile >> stringVal >> outputFreq;

    // convergence tolerance
    infile >> stringVal >> conv_tol;

    infile.close();

    cout << " Control parameters are successfully read " << endl;

    return;
}





void StabFEM::prepareInputData()
{
    printf("\n     StabFEM::prepareInputData()  .... STARTED ...\n");

    int ii, jj, kk, ee, nn, ind, n1, n2, dof;

    assert(ndim > 0 && ndim < 4);

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();

    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    // create elements and prepare element data
    elems = new ElementBase* [nElem];

    for(ee=0;ee<nElem;++ee)
    {
      if(ndim == 2)
      {
        //if(npElem == 3)
          //elems[ee] = new LagrangeElem2DNavierStokesTria3Node;
        //else if(npElem == 4)
          elems[ee] = new LagrangeElem2DNavierStokesQuad4Node;
      }
      else
      {
        //if(npElem == 4)
          //elems[ee] = new LagrangeElem3DNavierStokesTetra4Node;
        //else if(npElem == 8)
          //elems[ee] = new LagrangeElem3DNavierStokesHexa8Node;
      }

      elems[ee]->nodeNums = elemConn[ee];

      elems[ee]->prepareElemData(node_coords);
    }


    cout << " elements are created and prepared " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    soln.resize(nNode*ndof);
    soln.setZero();

    solnInit  = soln;
    solnPrev  = soln;
    solnPrev2 = soln;
    solnPrev3 = soln;
    solnPrev4 = soln;

    solnDot     = soln;
    solnDotPrev = soln;
    solnDotCur  = soln;
    solnApplied = soln;

    double  xx, yy, zz, fact;

    cout << " aaaaaaaaaaaaaaa " << endl;

    Kovasznay analy;

    solnApplied.setZero();
    for(ii=0; ii<nDBC; ++ii)
    {
        n1 = DirichletBCs[ii][0];
        n2 = DirichletBCs[ii][1];

        jj = n1*ndof+n2;

        //xx = node_coords[n1][0] ;
        //yy = node_coords[n1][1] ;

        //DirichletBCsVelo[ii][2] = analy.computeValue(n2, xx, yy);

        solnApplied(jj) = DirichletBCs[ii][2];
    }
    //printVector(solnApplied);

    printf("     StabFEM::prepareInputData()  .... FINISHED ...\n\n");

    return;
}




void StabFEM::assignBoundaryConditions(double timeCur, double dt, double timeFact)
{
    int ii, n1, n2, ind;
    for(ii=0; ii<nDBC; ++ii)
    {
        n1 = DirichletBCs[ii][0];
        n2 = DirichletBCs[ii][1];

        ind = n1*ndof+n2;

        solnApplied[ind] = DirichletBCs[ii][2] * timeFact - soln[ind];
        //cout << ii << '\t' << n1 << '\t' << n2 << '\t' << ind << '\t' << solnApplied[ind] << endl;

        //solnApplied[ind] = analy.computeValue(n2, node_coords[n1][0], node_coords[n1][1], 0.0, timeNow) - soln[ind];
    }
    //printVector(solnApplied);

    return;
}




void StabFEM::setInitialConditions()
{
    double  xx=0.0, yy=0.0, zz=0.0, fact;

    for(int ii=0; ii<nNode; ++ii)
    {
        xx = node_coords[ii][0];
        yy = node_coords[ii][1];
        zz = node_coords[ii][2];

        //veloPrev(ii*2) = 2.0*yy*(3.0-yy)/3.0;
        solnPrev(ii*ndim) = 1.0*yy;

        //solnPrev(ii*ndim) = 16.0*0.45*yy*zz*(0.41-yy)*(0.41-zz)/0.41/0.41/0.41/0.41;
    }
    soln = solnPrev;

    return;
}






void StabFEM::setTimeParam()
{
  //SolnData.setTimeParam();

  return;
}




void StabFEM::writeNodalData()
{
  return;
}





int StabFEM::initialise_pardiso()
{
  solnVec.resize(totalDOF);
  solnVec.setZero();

  perm.resize(totalDOF);

  phase = 11; error = 0;
  int  mtxType = 11, idum;

  char *tmp;


  mtxType =  1;  // real and structurally symmetric
  //mtxType = 11; // real and unsymmetric


  SOLVER = 0;       // sparse direct solver
  //SOLVER = 1;       // multi-recursive iterative solver
  MTYPE  = mtxType; // matrix type
  MAXFCT = 1;       // maximum number of factorisations of same sparsity pattern
                    //      to be kept in memory
  MNUM   = 1;       // which factorisation to use
  NRHS   = 1;       // number of right hand sides
  MSGLVL = 0;       // output message level (1 -> print statistical information)

  IPARM[0] = 0;     // PARADISO will set IPARM to default values
  //IPARM[0] = 1;     // user input values to IPARM
  //IPARM[1] = 2;

  tmp = getenv("OMP_NUM_THREADS");

  if (tmp != NULL) 
  {
    sscanf(tmp,"%d", &idum);
  }
  else printf("set environment variable OMP_NUM_THREADS!");


  cout << "OMP_NUM_THREADS = " << idum << endl;

  IPARM[2] = max(1,idum);  // number of processors (no default value available)


  pardisoinit_(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);

  if (error != 0)
  {
    if (error == -10) printf("no license file found.");
    if (error == -11) printf("license is expired.");
    if (error == -12) printf("wrong username or hostname.");
  }

  int  *c1, *c2, ii;

  csr.resize(totalDOF+1);
  col.resize(matK.nonZeros());

  c1 = matK.outerIndexPtr();
  c2 = matK.innerIndexPtr();


  for(ii=0;ii<=totalDOF;ii++)
    csr[ii] = c1[ii] + 1;

  for(ii=0;ii<matK.nonZeros();ii++)
    col[ii] = c2[ii] + 1;

  array = matK.valuePtr();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &totalDOF, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  if (error != 0)
  {
    cout << "PARDISO ERROR = " << error << "\n\n";
    printf("symbolic factorisation failed.");
  }

  //IPARM[4] = 0;  // user input permutation
  //IPARM[4] = 2;  // return the permutation
  IPARM[5] = 0; // do not overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  return 0;
}




int StabFEM::factoriseAndSolve_pardiso()
{
  phase = 23; error = 0;

  solnVec.setZero();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &totalDOF, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &solnVec[0], &error, DPARM);

  //cout << "error = " << error << endl;

//   printf("Peak memory [kB] phase 1       = %d \n", IPARM[14]);
//   printf("Permanent integer memory [kb]. = %d \n", IPARM[15]);
//   printf("Peak real memory [kB]          = %d \n", IPARM[16]);
//   printf("Number of nonzeros in LU.      = %d \n", IPARM[17]);
//   printf("Gflops for LU factorization.   = %d \n", IPARM[18]);

  return 0;
}
