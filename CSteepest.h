// ==================================================================//
// STEEPEST DESCENT METHOD
//
// IN THE MAIN PROGRAM, GIVE NUMBER OF VARIABLES, STARTING POINT,
// TOLERANCES, ITERATION LIMIT (you can experiment with STP & BET)
// IN SUBROUTINES GETFUN, GRADIENT, DEFINE YOUR FUNCTION AND ITS GRADIENT
// IN SUBROUTINE QUADFIT, GIVE THE VALUE OF THE TOLERANCE "EPS"
// ==================================================================//

#include <iostream.h>
#include <math.h>

const int MAX_NUM_VAR = 10;

class CSteepest
{

protected:
   int mNumVar;

   double mX [MAX_NUM_VAR];
   double mDF[MAX_NUM_VAR];
   double mD [MAX_NUM_VAR];
   double mX0[MAX_NUM_VAR];

   int    mIterLimit;
   double mTolGrad;
   double mEPSX;
   double mRelx;
   double mEPSF;
   double mRelf;


   // mNumVar = Number of variables
   // mStep   = Initial step to identify a 3-point pattern during line search
   // mBet    = Factor (< 1) by which mStep is reduced at each iteration
   // mEPSX   = A tolerance used in sub quadfit -- important!
   // -----------------------------------------------------------------------
   int    mNumVar;
   double mStep;
   double mBet;


public:
   CSteepest();


   double getFun(int& nfv);

   double getFun1(double al, int& nfv);

   double gradient();

}
