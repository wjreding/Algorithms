// ==================================================================//
// STEEPEST DESCENT METHOD
//
// IN THE MAIN PROGRAM, GIVE NUMBER OF VARIABLES, STARTING POINT,
// TOLERANCES, ITERATION LIMIT (you can experiment with STP & BET)
// IN SUBROUTINES GETFUN, GRADIENT, DEFINE YOUR FUNCTION AND ITS GRADIENT
// IN SUBROUTINE QUADFIT, GIVE THE VALUE OF THE TOLERANCE "EPS"
// ==================================================================//


#include <math.h>


// ==================================================================//
CSteepest::CSteepest()
{
   mIterLimit = 500;

   mTolGrad = 1.0e-6;
   mEPSX    = 1.0e-7;
   mRelx    = 1.0e-6;
   mEPSF    = 1.0e-6;
   mRelf    = 1.0e-6;

   mNumVar = 2;
   mStep   = 0.6;
   mBet    = 0.8;

   mX[0] = 1.0;
   mX[1] = 1.0;

}


// ==================================================================//
CSteepest::init(double x[], int size)
{
   mNumVar = size;

   for (int i = 0; i < size; i++)
   {
      mX[i] = x[i];
   }
}


// ==================================================================//
CSteepest::solveSteepest()
{
   int i;
   int nfv, iterm, ic;
   double f, x1, f1, xstar, fstar, fold, dfn;
   double fpr, ftol;


   for (i = 0; i < mNumVar; i++)
   {
      mX0[i] = mX[i];
   }

   nfv = 0;
   f   = getFun(nfv);

   x1    = 0;
   f1    = f;
   fstar = f;
   iterm = 0;
   ic    = 0;
   fold  = f;

   do 
   {
      iterm++;

      if (iterm > mIterLimit) 
      {
         cout << "*** iteration limit exceeded ***" << endl;
         break;
      }

      for (i = 0; i < mNumVar; i++)
      {
         mX[i] = mX0[i];
      }

      dfn = gradient();

      if (fabs(dfn) <= tolgrad)
      {
         break;
      }

      for (int i = 1; i <= n; i++) 
      {
         d[i] = -df[i] / dfn;
      }


      Quadfit(n, epsx, relx, epsf, relf, x0, x, d,x1, f1,xstar, fstar, stp, nfv);

      x1 = 0;
      f1 = fstar;

      for (i = 0; i < mNumVar; i++)
      {
         mX0[i] = mX0[i] + xstar * d[i];
      }


      mStepp  = mBbet * xstar;
      fpr     = fstar;
      ftol    = fabs(fstar) * mRelf + mEPSF;

      if (fabs(fpr - fold) < ftol) 
      {
         ic++;

         if (ic == 2)
            break;
      }
      else
      {
         ic = 0;
      }

      fold = fpr;

   } while (true);  // end do

}  // CSteepest::solveSteepest()


// ==================================================================//
//    ROSENBROCK'S FUNCTION
//    F = 100 * (X[1] * X[1] - X[2]) **2 + (1 - X[1]) **2
//    NONLINEAR SPRING, EXAMPLE 3.1 IN TEXT 
// ==================================================================//
double CSteepest::getFun(int& nfv)
{
   double ak1    = 1.0;
   double ak2    = 1.0;
   double force1 = 0.0;
   double force2 = 2.0;
   double delta1, delta2;
      
   nfv++;

   delta1 = sqrt( pow((mX[0] + 10.0), 2) + pow(mX[1] - 10.0), 2) ) - 10 * sqrt(2.0);

   delta2 = sqrt( pow((mX[0] - 10.0), 2) + pow((x[1] - 10.0), 2) ) - 10 * sqrt(2.0);

   f = (0.5* ak1 * pow(delta1, 2)) + (0.5 * ak2 * pow(delta2, 2)) - force1 * mX[0] - force2 * mX[1];

   return f;

}  // double CSteepest::getFun()


// ==================================================================//
double CSteepest::getFun1(double al, int& nfv)
{
   for (int i = 1; i <= mNumVar; i++)
   {
      mX[i] = mX0[i] + al * d[i];
   }

   f = getFun(nfv);

}  // double CSteepest::getFun1()


// ==================================================================//
// ROSENBROCK'S FUNCTION
// DF[1] =  400 * X[1] * (X[1] * X[1] - X[2]) - 2 * (1 - X[1])
// DF[2] = -200 * (X[1] * X[1] - X[2])
//
// NONLINEAR SPRING, EXAMPLE 3.1 IN TEXT 
// ==================================================================//
double CSteepest::gradient()
{
   double ak1    = 1.0;
   double ak2    = 1.0;
   double force1 = 0.0;
   double force2 = 2.0;

   double delta1, delta2;
   double dd1dx1, dd2dx1, dd1dx2, dd2dx2;

   delta1 = sqrt( pow((mX[0] + 10), 2) + pow((mX[1] - 10), 2) ) - 10 * sqrt(2.0);
   delta2 = sqrt( pow((mX[0] - 10), 2) + pow((mX[1] - 10), 2) ) - 10 * sqrt(2.0);

   dd1dx1 = (mX[0] + 10) / sqrt( pow((mX[0] + 10.0), 2) + pow((mX[1] - 10.0), 2) );
   dd2dx1 = (mX[0] - 10) / sqrt( pow((mX[0] - 10.0), 2) + pow((mX[1] - 10.0), 2) );
   dd1dx2 = (mX[1] - 10) / sqrt( pow((mX[0] + 10.0), 2) + pow((mX[1] - 10.0), 2) );
   dd2dx2 = (mX[1] - 10) / sqrt( pow((mX[0] - 10.0), 2) + pow((mX[1] - 10.0), 2) );

   mDF[0] = ak1 * delta1 * dd1dx1 + ak2 * delta2 * dd2dx1 - force1;
   mDF[1] = ak1 * delta1 * dd1dx2 + ak2 * delta2 * dd2dx2 - force2;

   // Compute norm of the gradient vector
   // -----------------------------------
   double dfn = 0.0;
   for (int i = 0; i < mNumVar; i++) 
   {
      dfn = dfn + mDF[i] * mDF[i];
   }

   dfn = sqrt(dfn);

   return dfn;

}  // double CSteepest::Gradient(double& dfn)


// ==================================================================//
	public static void QuadFit(int n, double epsx, double relx, double epsf, double relf,
	                           double [] x0, double [] xdv, double [] ddv, double x1, double f1,
	                           double xstar, double fstar, double stp, int nfv)
   {

		ThreePoint(n, x1, x2, x3, f1, f2, f3, stp, x0, xdv, ddv, nfv);
		
		a  = x1;
		b  = x3;
		fa = f1;
		fb = x3;
		sa = a;
		sb = b;

//    Need SA < SB
		if (b < a) {
			sa = b;
			sb = a;
		}

		x  = x2;
		fx = fb;
		w  = x;
		v  = w;
		e  = 0;
		fx = f2;
		fw = fx;
		fv = fw;

		do {
			sm  = .5 * (sa + sb);
			tol = relx * Math.Abs(x) + epsx;
			t2  = 2 * tol;

//       Convergence based on interval
			if ( Math.Abs(x - sm) <= (t2 - .5*(sb - sa)) ) exit;

			etmp = 0.0;
			p    = 0.0;
			q    = 0.0;

			if (Math.Abs(e) > tol) {
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2 * (q - r);
			
				if (q > 0.0) {
					p = -p;
				} else {
					q = -q;
				}
			
				etmp = e;
				e    = d;
			}  // if (Math.Abs(e) > tol)

			i1 = sign(1.0, q * (sa - x) - p);
			i2 = sign(1.0, q * (sb - x) - p);
			if (Math.Abs(p) >= Math.Abs(.5 * q * etmp) || i1 == i2) {
//				E is the length of the larger interval
				if (x < sm) {
					e = sb - x;
				} else {
					e = sa - x;
				}

				d = .381966 * e;

         } else {
				d = p / q;
				u = x + d;

				if (u - sa < t2 || sb - u < t2) {
					if (x < sm) {
						d = tol;
					} else { 
						d = -tol;
					}
				}
               
			}  // if (abs(p)

			if (Math.Abs(d) >= tol) {
				u = x + d;
			} else {
				if (d > 0.0) { 
					u = x + tol;
				} else { 
					u = x - tol;
				}
			}

			GetFun1(n, x0, xdv, ddv, u, fu, nfv);

//       Set A,B,X,U,V,W for the next iteration
			if (fu <= fx) {
				if (u < x) {
					sb = x;
					fb = fx;
				} else {
					sa = x;
					fa = fx;
				}

				v  = w;
				fv = fw;
				w  = x;
				fw = fx;
				x  = u;
				fx = fu;

			} else {
				if (u < x) {
					sa = u;
					fa = fu;
				} else {
					sb = u;
					fb = fu;
				}

				if ((fu <= fw) || (w == x)) {
					v  = w;
					fv = fw;
					w  = u;
					fw = fu;
				} else if ((fu <= fv) || (v == x) || (v == w)) {
					v  = u;
					fv = fu;
				}

			}  // if (fu <= fx)

//			Convergence based on function value
			tolf = relf * Math.Abs(fx) + epsf;
			if ( ((fa - fx) + (fb - fx)) < (2 * tolf) )
            break;

		}  while (true); // end do

		x2    = x;
		f2    = fx;
		xstar = x2;
		fstar = f2;

	}  // public void QuadFit


// ==================================================================//
   private static void threep() {
      double t;
      bool findBound;

      double tau = Math.Sqrt(1.25) - 0.5;

      nfv++;
      x1 = 0;
      f1 = GetFun(x1);

      nfv++;
      x2 = x1 + stp;
      f2 = GetFun(x2);


//      System.out.println();
//      System.out.println();
//      System.out.println("Number of Function Evaluations = " + nfv);
//      System.out.println("x1 = " + x1 + " x2 = " + x2 + " x3 = " + x3);
//      System.out.println("f1 = " + f1 + " f2 = " + f2 + " f3 = " + f3);


      if (f2 > f1) {
         t  = x1;
         x1 = x2;
         x2 = t;
         t  = f1;
         f1 = f2;
         f2 = t;
         stp = -stp;
      }

      findBound = true;
      while (findBound) {
         stp /= tau;
         x3 = x2 + stp;

         nfv++;
         f3 = getfun(x3);

         if (f3 > f2) {
            findBound = false;

         } else {
            f1 = f2;
            x1 = x2;
            f2 = f3;
            x2 = x3;
         }

      }  // while (findBound)

   }  // threep

}  // public class CSteepest 

