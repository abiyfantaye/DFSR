Info << "\nPerforming Cholesky decomposition for sampled frequencies" << nl << nl ; 

//The Cholesky decomposition is done using the following simplified form 
/*
  CPSD = diag(S)*C*diag(S)
       = diag(S)*Lc*Lc^T*diag(S)

  Where: 
        Lc: is the Cholesky factorization of the coherency matrix C
        S: is the auto-power spectrum at each point(size = n)
        diag(S) is a diagonal matrix created from S. 
  
  Thus, the Cholesky factorization of the CPSD matrix is computed as:

  L =  diag(S)*Lc

  Note: This approach prevents the Cholesky factorization from failing due to 
        non-positive semi-defnitness of the CPSD matrix caused by the inhomogeneity 
        in the coherency/covariance matrix 
*/

//Compute each entry of the coherency matrix.

vector coh(vector::zero);
vector s(vector::zero);

forAll(fInterp, f)
{  
  //Form and decompose the coherency matrix in 
  //the master processor.
  if(Pstream::master())
  {  
    forAll(faceCentres, i)
    {
      for (int j=0; j <= i; j++)
      {
        coh = Coh(fInterp[f], cohUav, C, faceCentres[j] - faceCentres[i]);
        Cu[j + i*(i+1)/2] = coh.x(); 
        Cv[j + i*(i+1)/2] = coh.y(); 
        Cw[j + i*(i+1)/2] = coh.z(); 
      }  
    }

    //Compute the Cholesky factorization using LAPACKE subroutine 
    //Uses a packed matrix storage.
    int infoU =  LAPACKE_dpptrf(LAPACK_ROW_MAJOR, 'L', n , Cu.data());
    int infoV =  LAPACKE_dpptrf(LAPACK_ROW_MAJOR, 'L', n , Cv.data());
    int infoW =  LAPACKE_dpptrf(LAPACK_ROW_MAJOR, 'L', n , Cw.data());

    //If decomposition fail, show error message and exit.
    if((infoU != 0)||(infoV != 0)||(infoW != 0))
    {
      FatalErrorInFunction
      << "ERROR: Cholesky decomposition failed!" << nl
      << exit(FatalError);
    }
  }

  Pstream::scatter(Cu);
  Pstream::scatter(Cv);
  Pstream::scatter(Cw);


  //Copy the decomposed matrix to H matrix in each processor. 
  forAll(facesIndx, i)
  {
    s = vonK(fInterp[f], mag(Uav[i]), I[i], L[i]);
    s.x() = ::sqrt(s.x());
    s.y() = ::sqrt(s.y());
    s.z() = ::sqrt(s.z());

    for (int j=0; j <= facesIndx[i]; j++)
    {
      H[i*n + j][f].x() = Cu[j + facesIndx[i]*(facesIndx[i]+1)/2]*s.x(); 
      H[i*n + j][f].y() = Cv[j + facesIndx[i]*(facesIndx[i]+1)/2]*s.y(); 
      H[i*n + j][f].z() = Cw[j + facesIndx[i]*(facesIndx[i]+1)/2]*s.z(); 
    }  
  }

  Info  << "Cholesky decomposition: " << 100*(f + 1.0)/nInterp << "\% completed" << ", f = " << fInterp[f] << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}

//Free up some space 
Cu.clear();
Cv.clear();
Cw.clear();
I.clear();
L.clear();