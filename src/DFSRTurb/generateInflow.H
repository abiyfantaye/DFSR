Info <<"\n\nSimulating velocity time history:  " << nl << nl ;

//Generated velocity at the current time. 
List<vectorField> allU(Nt, vectorField(nLoc, vector::zero));

//Number of frequencies in the FFT 
label M(2*nFreq);

//Arrays to hold inverse FFT variables.
fftw_complex *Bu, *Bv, *Bw, *Du, *Dv, *Dw;


Bu = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
Du = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

Bv = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
Dv = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

Bw = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
Dw = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

//FFTW plans 
fftw_plan pu;
fftw_plan pv;
fftw_plan pw;

scalar sqrt2df = ::sqrt(2.0*df);
vector Hjm;
scalar f,temp, cs, sn;
label q;

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

                    Bu[l][0] = Hjm.x()*cosPhi[l*n + m].x();
                    Bv[l][0] = Hjm.y()*cosPhi[l*n + m].y();    
                    Bw[l][0] = Hjm.z()*cosPhi[l*n + m].z();
                    

                    Bu[l][1] = Hjm.x()*sinPhi[l*n + m].x();
                    Bv[l][1] = Hjm.y()*sinPhi[l*n + m].y();
                    Bw[l][1] = Hjm.z()*sinPhi[l*n + m].z();
               }
               else
               {
                    Bu[l][0] = 0.0; Bu[l][1] = 0.0;
                    Bv[l][0] = 0.0; Bv[l][1] = 0.0;
                    Bw[l][0] = 0.0; Bw[l][1] = 0.0;
               }
          }

          pu = fftw_plan_dft_1d(M, Bu, Du, FFTW_BACKWARD, FFTW_ESTIMATE);
          pv = fftw_plan_dft_1d(M, Bv, Dv, FFTW_BACKWARD, FFTW_ESTIMATE);
          pw = fftw_plan_dft_1d(M, Bw, Dw, FFTW_BACKWARD, FFTW_ESTIMATE);

          fftw_execute(pu);
          fftw_execute(pv);
          fftw_execute(pw);

          fftw_destroy_plan(pu);
          fftw_destroy_plan(pv);
          fftw_destroy_plan(pw);

          f = dt*twopi*(m + 1.0)*df/n;

          for (label t = 0; t < Nt; t++)
          {
               temp = t*f;
               cs = ::cos(temp);
               sn = ::sin(temp);
               q = t%M;

               allU[t][i].x() += (Du[q][0]*cs - Du[q][1]*sn); 
               allU[t][i].y() += (Dv[q][0]*cs - Dv[q][1]*sn); 
               allU[t][i].z() += (Dw[q][0]*cs - Dw[q][1]*sn); 
          }  
     } 

     Info << "Generation completed for "<< i*Pstream::nProcs() << " points out of " << n << nl
          << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << nl << endl;
    
}

//Clear variables that are not needed anymore
cosPhi.clear();
sinPhi.clear();
H.clear();

//Free FFTW space 
fftw_free(Bu); fftw_free(Du); 
fftw_free(Bv); fftw_free(Dv); 
fftw_free(Bw); fftw_free(Dw);

//Add meant to the fluctuation 
for(label t = 0; t < Nt; t++)
{    
  allU[t] += Uav; 
}