//Write only from master processor 
if(Pstream::master())
{    
     if(writeInflow)
     {
          fileName timePath = inflowPath/runTime.timeName();
          mkDir(timePath);

          OFstream uOS(timePath/"U");
          uOS << Ut;
     }

     if(writeSamples)
     {                
          sampleUxOS.open(samplePath/"Ux", std::ios_base::app);
          sampleUyOS.open(samplePath/"Uy", std::ios_base::app);
          sampleUzOS.open(samplePath/"Uz", std::ios_base::app);   

          sampleUxOS <<  runTime.timeName() << "\t";
          sampleUyOS <<  runTime.timeName() << "\t";
          sampleUzOS <<  runTime.timeName() << "\t";

          forAll(samplingPoints, i)
          {
               sampleUxOS << Ut[samplingPointsIndx[i]].x() << "\t";
               sampleUyOS << Ut[samplingPointsIndx[i]].y() << "\t";
               sampleUzOS << Ut[samplingPointsIndx[i]].z() << "\t";
          }

          sampleUxOS << nl;
          sampleUyOS << nl;
          sampleUzOS << nl;  

          sampleUxOS.close();
          sampleUyOS.close();
          sampleUzOS.close();     

          if(runTime.timeIndex()%(Nt/nVtkSamples) == 0)
          {                    
               Info <<"Writting the velocity field sample to a VTK file" << nl;

               autoPtr<surfaceWriter> formater = surfaceWriter::New("vtk", IOstream::ASCII);

               formater->write
               (
                    samplePath/"vtkSlices/",   // <case>/surface/TIME
                    runTime.timeName(),                      // name of surface
                    mergedPoints,
                    mergedFaces,
                    "Ut",             // name of field
                    Ut,
                    false               //point(true) or cell(false) data
               );
          }
     }
}