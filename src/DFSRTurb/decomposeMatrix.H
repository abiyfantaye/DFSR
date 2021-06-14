Info << "\nPerforming Cholesky decomposition for sampled frequencies" << nl << nl ; 

//The Cholesky decomposition is done using the following simplified form 
/*
  CPSD = diag(S)*Coh*diag(S)
       = diag(S)*Lc*Lc^T*diag(S)

  Where: 
        Lc: is the Cholesky factorization of the coherency matrix Coh
        S: is the auto-power spectrum at each point(size = n)
        diag(S) is a diagonal matrix created from S. 
  
  Thus, the Cholesky factorization of the CPSD matrix is computed as:

  L =  diag(S)*Lc

  Note: This approach prevents the Cholesky factorization from failing due to 
        non-positive semi-defnitness of the CPSD matrix caused by the inhomogeneity 
        in the coherency/covariance matrix 
*/

//Compute each entry of the coherency matrix.

forAll(fInterp, f)
{  
  //Form and decompose the coherency matrix in 
  //the master processor.

  for(label cmpt=0; cmpt < 3; cmpt++)
  {
    if(Pstream::master())
    {  
      forAll(faceCentres, i)
      {
        for (int j=0; j <= i; j++)
        {
          Lc[j + i*(i+1)/2] = Coh(fInterp[f], cohUav, C, faceCentres[j] - faceCentres[i], cmpt); 
        }  
      }

      //Compute the Cholesky factorization using LAPACKE subroutine 
      //Uses a packed matrix storage.
      int info =  LAPACKE_dpptrf(LAPACK_ROW_MAJOR, 'L', n , Lc.data());

      //If decomposition fail, show error message and exit.
      if(info)
      {
        FatalErrorInFunction
        << "ERROR: Cholesky decomposition failed for component: " << cmpt << nl
        << exit(FatalError);
      }
    }

    Pstream::scatter(Lc);

    //Copy the decomposed matrix to H matrix in each processor. 
    forAll(facesIndx, i)
    {
      for (int j=0; j <= facesIndx[i]; j++)
      {
        H[i*n + j][f][cmpt] = Lc[j + facesIndx[i]*(facesIndx[i]+1)/2]*::sqrt(vonK(fInterp[f], mag(Uav[i]), I[i][cmpt], L[i][cmpt], cmpt)); 
      }  
    }

  }

  Info  << "Cholesky decomposition: " << 100*(f + 1.0)/nInterp << "\% completed" << ", f = " << fInterp[f] << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}

//Free up some space 
Lc.clear();
I.clear();
L.clear();