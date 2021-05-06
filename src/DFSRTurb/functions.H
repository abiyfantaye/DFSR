//- Returns the von Karman spectrum for i-th component given the
//  mean velocity(Uav), turbulence intensity(I) and length scales(L)
inline scalar vonK
(
  const label cmpt,
  const scalar f,
  const scalar Uav,
  const scalar I,
  const scalar L
)
{
  if(cmpt==0){
    return 4.0*::pow(I*Uav, 2.0)*(L/Uav)/::pow(1.0 + 70.8*::pow(f*L/Uav, 2.0), 5.0/6.0);
  }
  else if(cmpt==1){
    return 4.0*::pow(I*Uav, 2.0)*(L/Uav)*(1.0 + 188.4*::pow(2.0*f*L/Uav, 2.0))/::pow(1.0 + 70.8*::pow(2.0*f*L/Uav, 2.0), 11.0/6.0);
  }
  else if(cmpt==2){
    return 4.0*::pow(I*Uav, 2.0)*(L/Uav)*(1.0 + 188.4*::pow(2.0*f*L/Uav, 2.0))/::pow(1.0 + 70.8*::pow(2.0*f*L/Uav, 2.0), 11.0/6.0);
  }

  FatalErrorInFunction
  << "Component i = " << cmpt << " is unknown \n"
  << exit(FatalError);
  return 0;
}


//- Returns the coherency function for i-th component given freqquency,
//  coherency decay coefficient(C), mean velocity(Uav), and point pair(p1, p2).
inline scalar Coh
(  
    const label cmpt, 
    const scalar f,     
    const scalar Uav, 
    const tensor& C, 
    const vector& p1, 
    const vector& p2
)
{
  if(cmpt==0){
    return ::exp(-f*::sqrt(C.xx()*C.xx()*::pow(p2.x() - p1.x(), 2.0) + C.xy()*C.xy()*::pow(p2.y() - p1.y(), 2.0) + C.xz()*C.xz()*::pow(p2.z() - p1.z(), 2.0))/Uav);  
  }  
  else if(cmpt==1){
    return ::exp(-f*::sqrt(C.yx()*C.yx()*::pow(p2.x() - p1.x(), 2.0) + C.yy()*C.yy()*::pow(p2.y() - p1.y(), 2.0) + C.yz()*C.yz()*::pow(p2.z() - p1.z(), 2.0))/Uav); 
  }
  else if(cmpt==2){
    return ::exp(-f*::sqrt(C.zx()*C.zx()*::pow(p2.x() - p1.x(), 2.0) + C.zy()*C.zy()*::pow(p2.y() - p1.y(), 2.0) + C.zz()*C.zz()*::pow(p2.z() - p1.z(), 2.0))/Uav); 
  }
  
  FatalErrorInFunction
  << "Component i = " << cmpt << " is unknown \n"
  << exit(FatalError);
  return 0;
}

void printFullMatrix(label n, scalarField& matrix)
{
  printf("\n \n");
  for (int i=0; i<n; i++)
  {
    for(int j=0; j<n; j++) 
    {
      printf("%.4f\t",matrix[i*n + j]);
    }
    printf("\n");
  }
}

//Print the lower triangular part of a packed matrix(double)
void printPackedMatrix(label n, scalarField& matrix)
{
  printf("\n \n");
  for (int i=0; i<n; i++)
  {
    for(int j=0; j<=i; j++) 
    {
      printf("%.4f\t", matrix[j + i*(i+1)/2]);
    }
    printf("\n");
  }
}

scalar interpolateSpline
(
    scalar x,
    scalarField& xOld,
    scalarField& yOld
)
{
    label n = xOld.size();

    // early exit if out of bounds or only one value
    if (n == 1 || x < xOld[0])
    {
        return yOld[0];
    }
    if (x > xOld[n - 1])
    {
        return yOld[n - 1];
    }

    // linear interpolation if only two values
    if (n == 2)
    {
        return (x - xOld[0])/(xOld[1] - xOld[0])*(yOld[1] - yOld[0]) + yOld[0];
    }

    // find bounding knots
    label hi = 0;
    while (hi < n && xOld[hi] < x)
    {
        hi++;
    }

    label lo = hi - 1;


    //Info << "x = " << x << nl

    //Info << "lo = " << xOld[lo] << nl;
    //Info << "hi = " << xOld[hi] << nl;

    const scalar y1 = yOld[lo];
    const scalar y2 = yOld[hi];

    scalar y0;
    if (lo == 0)
    {
        y0 = 2*y1 - y2;
    }
    else
    {
        y0 = yOld[lo - 1];
    }

    scalar y3;
    if (hi + 1 == n)
    {
        y3 = 2*y2 - y1;
    }
    else
    {
        y3 = yOld[hi + 1];
    }

    // weighting
    scalar mu = (x - xOld[lo])/(xOld[hi] - xOld[lo]);

    // interpolate
    return
        0.5
       *(
            2*y1
          + mu
           *(
               -y0 + y2
              + mu*((2*y0 - 5*y1 + 4*y2 - y3) + mu*(-y0 + 3*y1 - 3*y2 + y3))
            )
        );
}
