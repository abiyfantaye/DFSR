//Create streams and paths 
fileName inflowPath = runTime.constant()/"boundaryData"/patchName;
fileName samplePath = runTime.constant()/"boundaryData"/"windProfile"/"sampledData";

vectorField samplingPoints;
labelList samplingPointsIndx;

ofstream sampleUxOS;
ofstream sampleUyOS;
ofstream sampleUzOS;


const label nVtkSamples = 4;

//Create the directory where the inflow 
//will be stored. 
if(Pstream::master())
{
     if(writeInflow)
     {
          mkDir(inflowPath);

          OFstream pointsOs(inflowPath/"points");
          pointsOs << faceCentres; 
     }

     if(writeSamples)
     {
          mkDir(samplePath);
          mkDir(samplePath/"vtkSlices");
          
          IFstream sampleIS(samplePath/"samplingPoints");
          sampleIS >> samplingPoints; 

          samplingPointsIndx = labelList(samplingPoints.size(), -1);
          scalar dist = VGREAT;

          //Find the index of the nearest points to the samplingPoints
          for(label i = 0; i < samplingPoints.size(); i++)
          {
               for(int j = 0; j < n ; j++)
               {
                    if(mag(faceCentres[j] - samplingPoints[i]) < dist)
                    {
                         dist = mag(faceCentres[j] - samplingPoints[i]);
                         samplingPointsIndx[i] = j;
                    }
               }
               dist = VGREAT;
          }
          
          sampleUxOS.open(samplePath/"Ux");
          sampleUyOS.open(samplePath/"Uy");
          sampleUzOS.open(samplePath/"Uz");

          // //Write the nearest points 
          if(sampleUxOS.good() && sampleUyOS.good() && sampleUzOS.good())
          {
            sampleUxOS << "# Nearest points (x, y, z)" << nl;
            sampleUyOS << "# Nearest points (x, y, z)" << nl;
            sampleUzOS << "# Nearest points (x, y, z)" << nl;          
            
            forAll(samplingPoints, i)
            {
                sampleUxOS << "# Probe " << i << " ("   << faceCentres[samplingPointsIndx[i]].x() << ", " 
                                                        << faceCentres[samplingPointsIndx[i]].y() << ", " 
                                                        << faceCentres[samplingPointsIndx[i]].z() << ")" << nl;

                sampleUyOS << "# Probe " << i << " ("   << faceCentres[samplingPointsIndx[i]].x() << ", " 
                                                        << faceCentres[samplingPointsIndx[i]].y() << ", " 
                                                        << faceCentres[samplingPointsIndx[i]].z() << ")" << nl;

                sampleUzOS << "# Probe " << i << " ("   << faceCentres[samplingPointsIndx[i]].x() << ", " 
                                                        << faceCentres[samplingPointsIndx[i]].y() << ", " 
                                                        << faceCentres[samplingPointsIndx[i]].z() << ")" << nl;
            }

            sampleUxOS << "#Time " << "\t";
            sampleUyOS << "#Time " << "\t";
            sampleUzOS << "#Time " << "\t";

            forAll(samplingPoints, i)
            {
                sampleUxOS << i << "\t";
                sampleUyOS << i << "\t";
                sampleUzOS << i << "\t";
            }

            sampleUxOS << nl;
            sampleUyOS << nl;
            sampleUzOS << nl;
          }

          sampleUxOS.close();
          sampleUyOS.close();
          sampleUzOS.close();
     }

}
