//- Returns the von Karman spectrum for component: cmpt 
//  mean velocity(Uav), turbulence intensity(I) and length scales(L)
inline scalar vonK
(
  const scalar f,
  const scalar Uav,
  const scalar I,
  const scalar L,
  const label cmpt

)
{
    if(cmpt==0)
    {
        return 4*::pow(I*Uav, 2)*(L/Uav)/::pow(1 + 70.8*::pow(f*L/Uav, 2), 5.0/6.0);
    }
    else
    {
        return 4*::pow(I*Uav, 2)*(L/Uav)*(1 + 188.4*::pow(2*f*L/Uav, 2))/::pow(1.0 + 70.8*::pow(2*f*L/Uav, 2), 11.0/6.0);
    }
}


//- Returns the coherency vector for a given freqquency,
//  coherency decay coefficient(C), mean velocity(Uav), and deference between points(p1, p2).
inline scalar Coh
(  
  const scalar f,     
  const scalar Uav, 
  const tensor& C, 
  const vector& dp,
  const label cmpt 
)
{
  return ::exp(-f*::mag(cmptMultiply(C.vectorComponent(cmpt), dp))/Uav);  
}

//Interplates a given vectorField with cubic splines 
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
