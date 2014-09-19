double average_rgyr(double *a, int start, int end) //function declaration

{ // start is the first index to use, end the end of the array, max index

double m; //declrations local to function
double s;

int n;

s=0.0;
for (n = start; n <= end; n++)
  {
     s += a[n];
  }
m = (s / n);

return (m);
}


