Info << "\nPerforming Cholesky decomposition for sampled frequencies" << nl << nl ; 

//The CPSD matrix is symmetric and only the lower triangular part is stored.
//The matrix is stored in row-major packed Lapack format with n*(n+1)/2 entries 
//including the diagonal elements
scalarField cpsd(n*(n+1)/2, 0.0);

//Initialize the spectrum at any particular frequency
scalarField S(n, 0.0);


label cmpt = 2;

//The Cholesky decomposition is done using the following simplified form 
/*
  CPSD = diag(S)*C*diag(S)
  CPSD = diag(S)*Lc*Lc^T*diag(S)

  Where: 
        Lc: is the Cholesky factorization of the coherency matrix C
        S: is the auto-power spectrum at each point(size = n)
        diag(S) is a diagonal matrix created from S. 
  
  Thus, the Cholesky factorization of the CPSD matrix is computed as:

  L =  diag(S)*Lc

  Note: This approach prevents the Cholesky factorization from failing due to 
        non-positive semi-defnitness of the CPSD matrix caused by the inhomogeneity 
        in the coherency matrix(covariance matrix) 
*/

//Compute each entry of the coherency matrix.
forAll(fInterp, f)
{  
  forAll(points, i)
  {
    S[i] = vonK(cmpt, fInterp[f], Uav[i], I[i][cmpt], L[i][cmpt]);

    for (int j=0; j <= i; j++)
    {
      cpsd[j + i*(i+1)/2] = Coh(cmpt, fInterp[f], cohUav, C, points[i], points[j]); 
    }  
  }

  //Compute the Cholesky factorization using LAPACKE subroutine 
  //Uses a packed matrix storage.
  int info =  LAPACKE_dpptrf(LAPACK_ROW_MAJOR, 'L', n , cpsd.data());

  //If decomposition fail, show error message and exit.
  if((info != 0))
  {
    FatalErrorInFunction
    << "ERROR: Cholesky decomposition failed for component = " << cmpt << nl
    << exit(FatalError);
  }

  S = sqrt(S);

  //Copy the decomposed matrix to H matrix 
  forAll(points, i)
  {
    for (int j=0; j <= i; j++)
    {
      H[j + i*(i+1)/2][f] = cpsd[j + i*(i+1)/2]*S[i]; 
    }  
  }

  Info  << "Cholesky decomposition" << 100*(f + 1.0)/nInterp << "\% completed" << ", f = " << fInterp[f] << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
}

cpsd.clear();
S.clear();
 
