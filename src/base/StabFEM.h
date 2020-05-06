
#ifndef incl_StabFEM_CLASS_h
#define incl_StabFEM_CLASS_h
#define EIGEN_SUPERLU_SUPPORT

#include <string.h>
#include <vector>
#include <fstream>
#include "ElementBase.h"
#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>


using std::vector;
using std::cout;
using std::string;
using std::ofstream;

class ElementBase;


#define ELEMENT_TYPE_NAMES_2D {\
                            "LagrangeElem2DNavierStokesTria3Node", \
                            "LagrangeElem2DNavierStokesQuad4Node", NULL};

#define ELEMENT_TYPE_NAMES_3D {\
                            "LagrangeElem3DNavierStokesTetra4Node", \
                            "LagrangeElem3DNavierStokesHexa8Node", NULL};

#define ELEMENT_TYPE_NAMES_FACE {\
                                "BernsteinElem2DEdge3Node", \
                                "BernsteinElem3DFaceTria6Node", \
                                "BernsteinElem3DFaceQuad9Node", NULL};

enum TimeStepType
{
  EXPLICIT = 0, SEMI_IMPLICIT = 1, FULL_IMPLICIT = 2
};


class StabFEM
{
    //private:

    public:

        int  ndim, ndof, npElem, nElem, totalDOF;
        int  nNode, nDBC, nFBC, nOutputFaceLoads, fileCount;
        int  stepsMax, outputFreq, tis;
        int  AlgoType;

        double  conv_tol, rhoInf, timeFinal, dt;
        double  elemData[50], timeData[50];

        vector<vector<double> >  node_coords;               //!< coordinates of the nodes (or control points)
        vector<vector<int> >     outputEdges;               //!< data for computing drag/lift forces
        vector<vector<int> >     elemConn;                  //!< element-node connectivity array

        vector<int>  assyForSoln, OutputNodes;

        vector<vector<double> >  DirichletBCs;          //!< Dirichlet BCs
        vector<vector<double> >  NeumannBCs;                //!< Neumann BCs
        vector<vector<double> >  InitialConds;              //!< Initial conditions
        vector<vector<double> >  OutputData;                //!< data for output
        vector<vector<double> >  nodeForcesData;
        vector<vector<double> >  ElemFaceLoadData;


        ElementBase  **elems;
        ElementBase  **elemsFaces;

        VectorXd  solnInit, ForceVectorExternal;
        VectorXd  totalForce, totalMoment, centroid;
        VectorXd  solnDot, solnDotPrev, solnDotCur;
        VectorXd  soln, solnCur, solnPrev, solnPrev2, solnPrev3, solnPrev4;
        VectorXd  solnApplied, rhsVec, solnVec;

        string  infilename;
        ofstream  fout_convdata;

        SparseMatrixXd  matK;

        //SimplicialLDLT<SparseMatrix<double> > solver;
        //SparseLU<SparseMatrixXd > solver;
        //SuperLU<SparseMatrixXd > solver;
        BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

        int   PT[64], IPARM[64];
        int   phase, error, SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL;

        double DPARM[64], ddum, *array;

        vector<int>  csr, col, perm;

    public:

        StabFEM();

        ~StabFEM();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  assignBoundaryConditions(double timeCur, double dt, double timeFact);

        void  prepareInputData();

        void  readInputData(string& fname);

        void  readControlParameters(string& fname);

        void  writeNodalData();

        void  writeReadResult(int, string&);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int   setSolver(int, int *parm = NULL, bool cIO = false);

        int   prepareMatrixPattern();

        int   setSolverDataForFullyImplicit();

        void  setTimeParam();

        void  timeUpdate();

        void  addExternalForces(double loadFact);

        void  computeElementErrors(int);

        void  setInitialConditions();

        int   solveFullyImplicit();

        int   initialise_pardiso();

        int   factoriseAndSolve_pardiso();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  postProcess();
};






#endif






