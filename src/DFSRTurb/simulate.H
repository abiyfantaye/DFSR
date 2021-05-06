Info <<"\n\nSimulating velocity time history:  " << nl << nl ;

label M(2*nFreq);

fftw_complex *B, *D;

B = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
D = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

fftw_plan plan;
scalar pi = Foam::constant::mathematical::pi;

//Simulate the time series using FFT. 
forAll(points, j)
{
     for (int m = 0; m <= j; m++)
     {
          for(int l = 0; l < M; l++)
          {
               if(l < nFreq)
               { 
                    scalar Hjm = interpolateSpline((l + (m + 1.0)/n)*df, fInterp, H[j*n + m]);
                    B[l][0] = Hjm*::sqrt(2.0*df)*::cos(Phi[l*n + m][cmpt]); 
                    B[l][1] = Hjm*::sqrt(2.0*df)*::sin(Phi[l*n + m][cmpt]);
               }
               else
               {
                    B[l][0] = 0.0; 
                    B[l][1] = 0.0;
               }
          }

          plan = fftw_plan_dft_1d(M, B, D, FFTW_BACKWARD, FFTW_ESTIMATE);

          fftw_execute(plan);
          fftw_destroy_plan(plan);

          for (int p = 0; p < Nt; p++)
          {
               U[p][j][cmpt] += ( D[p%M][0]*::cos(p*dt*2*pi*(m + 1.0)*df/n) 
                                - D[p%M][1]*::sin(p*dt*2*pi*(m + 1.0)*df/n)); 
          }
     }
}

fftw_free(B); 
fftw_free(D); 

//Pout <<"Processor "<< Pstream::myProcNo() << " finished simulation!" << nl << nl; 
