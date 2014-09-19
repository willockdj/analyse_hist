#include <math.h>
double std_rgyr(double *a, int start, int end)

{ // i is number of elements in array
         double s = 0;
         double m =0;
         double d = 0;
         double si =0;  //declrations local to function
         int i, n;

         for (n = start; n <= end; n++)

         {
               s  = s + a[n];
         }
          m = (s / n);
          for (n = start; n <= end; n++)
           {
               si = si + pow(a[n]-m,2);
           }
          d = sqrt(si/(n-1));
         return (d);
}


