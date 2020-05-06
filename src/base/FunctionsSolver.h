
#ifndef incl_FunctionsSolver_h
#define incl_FunctionsSolver_h


extern "C"
{
  void pardisoinit_(void*, int*, int*, int*, double*, int*);

  void pardiso_(void*, int*, int*, int*, int*, int*, double*,
                int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);
}


#endif

