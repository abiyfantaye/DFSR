//Fill the double indexed phase angles
if(Pstream::master())
{      
   //Pick a random phase angels from a uniform distribution 
   //in the interval [0.0, 2*pi]
   for(int m = 0; m < n ; m++)
   {
      for(int l = 0; l < nFreq ; l++)
      {
         sinPhi[m*nFreq + l]  =  twopi*rand.sample01<scalar>(); 
                  
         cosPhi[m*nFreq + l]  =  ::cos(sinPhi[m*nFreq + l]); 
         sinPhi[m*nFreq + l]  =  ::sin(sinPhi[m*nFreq + l]); 
      }
   }
}

Pstream::scatter(sinPhi);
Pstream::scatter(cosPhi);

//Simulate the velocity at each point using FFT
for (label i = 0; i < nLoc; i++)
{
     label j = facesIndx[i];

     for (label m = 0; m <= j; m++)
     {    
          for(label l = 0; l < M; l++)
          {
               if(l < nFreq)
               { 
                    Hjm = sqrt2df*interpolateSpline((l + (m + 1.0)/n)*df, fInterp, H[i*n + m]);

                    B[l][0] = Hjm*cosPhi[l*n + m];
                    B[l][1] = Hjm*sinPhi[l*n + m];

               }
               else
               {
                    B[l][0] = 0.0; 
                    B[l][1] = 0.0;
               }
          }

          pln = fftw_plan_dft_1d(M, B, D, FFTW_BACKWARD, FFTW_ESTIMATE);

          fftw_execute(pln);

          fftw_destroy_plan(pln);

          freq = dt*twopi*(m + 1.0)*df/n;

          for (label t = 0; t < Nt; t++)
          {
               allU[t][i][cmpt] += (D[t%M][0]*::cos(t*freq) - D[t%M][1]*::sin(t*freq)); 
          }  
     } 

     Info << "Generating for " << cmptName[cmpt] << " : "<< (i + 1)*Pstream::nProcs() << " points out of " << n << nl
          << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << nl << endl;
    
}

if(cmpt == nCmpt-1)
{
     //Clear variables that are not needed anymore
     cosPhi.clear();
     sinPhi.clear();
     H.clear();

     //Free FFTW space 
     fftw_free(B); 
     fftw_free(D); 

     //Add meant to the fluctuation 
     for(label t = 0; t < Nt; t++)
     {    
          allU[t] += Uav; 
     }
}