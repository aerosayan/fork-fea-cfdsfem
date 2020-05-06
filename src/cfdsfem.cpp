/*
! Program for 2D explicit finite element analysis of incompressible Navier-Stokes
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 17-May-2018
! Place : Swansea, UK
!
!
*/


#include "headersVTK.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "StabFEM.h"

using namespace std;



int main(int argc, char* argv[])
{
    //Set the input file name
    //The file name is specified from the command line
    if(argc == 0)
    {
        cerr << " Error in input data " << endl;
        cerr << " You must enter name of input file" << endl;
        cerr << " Aborting..." << endl;
    }

    string  meshfile = argv[1];
    string  controlfile = argv[2];


    StabFEM  stabfem;

    stabfem.readInputData(meshfile);

    stabfem.readControlParameters(controlfile);

    stabfem.prepareInputData();

    stabfem.solveFullyImplicit();

    //if( !infileOutput.empty() )

    //cout << " cccccccccccccc " << endl;

    cout << " Program is successful \n " << endl;

    return 0;
}

