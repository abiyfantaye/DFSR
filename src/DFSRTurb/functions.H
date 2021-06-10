//- Returns the von Karman spectrum all components given the
//  mean velocity(Uav), turbulence intensity(I) and length scales(L)
inline vector vonK
(
  const scalar f,
  const scalar Uav,
  const vector I,
  const vector L
)
{
  vector s(vector::zero);

  s.x() = 4.0*::pow(I.x()*Uav, 2.0)*(L.x()/Uav)/::pow(1.0 + 70.8*::pow(f*L.x()/Uav, 2.0), 5.0/6.0);
  s.y() = 4.0*::pow(I.y()*Uav, 2.0)*(L.y()/Uav)*(1.0 + 188.4*::pow(2.0*f*L.y()/Uav, 2.0))/::pow(1.0 + 70.8*::pow(2.0*f*L.y()/Uav, 2.0), 11.0/6.0);
  s.z() = 4.0*::pow(I.z()*Uav, 2.0)*(L.z()/Uav)*(1.0 + 188.4*::pow(2.0*f*L.z()/Uav, 2.0))/::pow(1.0 + 70.8*::pow(2.0*f*L.z()/Uav, 2.0), 11.0/6.0);

  return s;
}


//- Returns the coherency vector for a given freqquency,
//  coherency decay coefficient(C), mean velocity(Uav), and deference between points(p1, p2).
inline vector Coh
(  
  const scalar f,     
  const scalar Uav, 
  const tensor& C, 
  const vector& dp
)
{
  vector coh(vector::zero);
  coh.x() = ::exp(-f*mag(vector(C.xx()*dp.x(), C.xy()*dp.y(), C.xz()*dp.z()))/Uav);  
  coh.y() = ::exp(-f*mag(vector(C.yx()*dp.x(), C.yy()*dp.y(), C.yz()*dp.z()))/Uav); 
  coh.z() = ::exp(-f*mag(vector(C.zx()*dp.x(), C.zy()*dp.y(), C.zz()*dp.z()))/Uav); 
  
  return coh;
}

//Interplates a given vectorField with cubic splines 
vector interpolateSpline
(
    scalar x,
    scalarField& xOld,
    vectorField& yOld
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

    const vector& y1 = yOld[lo];
    const vector& y2 = yOld[hi];

    vector y0;
    if (lo == 0)
    {
        y0 = 2*y1 - y2;
    }
    else
    {
        y0 = yOld[lo - 1];
    }

    vector y3;
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