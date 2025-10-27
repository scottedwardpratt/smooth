#include "msu_smooth/master.h"
using namespace std;
using namespace NBandSmooth;

void CSmoothMaster::TestAtTrainingPts(){
	unsigned int iY;
	for(iY=0;iY<observableinfo->NObservables;iY++){
		TestAtTrainingPts(iY);
	}

}

void CSmoothMaster::TestAtTrainingPts(unsigned int iY){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain;
	double Y,SigmaY_emulator;
	CLog::Info("----------- "+observableinfo->observable_name[iY]+" -------------\n");
	CLog::Info("---- Lambda="+to_string(emulator[iY]->LAMBDA)+", SigmaA="+to_string(emulator[iY]->SigmaA)+"\n");
   CLog::Info("- itrain --- Y_full     Y_emulator    Sigma_emulator\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info(" "+to_string(itrain)+" ");
		GetY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e  +/- %12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::TestAtTrainingPts(string obsname){
	char pchars[CLog::CHARLENGTH];
	unsigned int itrain,iY;
	double Y,SigmaY_emulator;
	iY=observableinfo->GetIPosition(obsname);
	CLog::Info("--- TESTING AT TRAINING POINTS ----\n");
	for(itrain=0;itrain<traininginfo->NTrainingPts;itrain++){
		CLog::Info("------ itrain="+to_string(itrain)+" --------\n");
		GetY(iY,traininginfo->modelpars[itrain],Y,SigmaY_emulator);
		snprintf(pchars,CLog::CHARLENGTH,
		"Y[%u]=%10.3e =? %10.3e,    SigmaY=%12.5e\n",iY,traininginfo->YTrain[iY][itrain],Y,SigmaY_emulator);
		CLog::Info(pchars);
	}
}

void CSmoothMaster::TestVsFullModel(){
   char cdummy[200];
	unsigned int iY,ntest=0,itest;
	unsigned int NObs=observableinfo->NObservables,NYY;
	vector<double> Y,SigmaY_emulator;
	vector<double> testtheta;
   vector<FILE *> fptr_YvsY;
   FILE *fptr_in;
	string filename;
	fitpercentage=0.0;
   vector<int> nfit;
   vector<double> averagesigma2;
   string(command);
   double realY;
   NYY=NObs;
   // YOu can't open too many file pointers!
   if(NYY>=100)
      NYY=100;
	
	int ifit;
	int nfitpercent[100]={0};
   ReadTestingInfo();
   ntest=CSmoothEmulator::NTestingPts;
   nfit.resize(NObs,0);
   averagesigma2.resize(NObs,0.0);
   Y.resize(NObs);
   SigmaY_emulator.resize(NObs);
   fptr_YvsY.resize(NObs);
   command="mkdir -p smooth_data/output_stuff/fullmodel_testdata";
   system(command.c_str());
   for(iY=0;iY<NYY;iY++){
      filename="smooth_data/output_stuff/fullmodel_testdata/YvsY_"+observableinfo->observable_name[iY]+".txt";
      fptr_YvsY[iY]=fopen(filename.c_str(),"w");
   }
   
   for(itest=0;itest<ntest;itest++){
      GetAllY(testinginfo->modelpars[itest],Y,SigmaY_emulator); // emulated values
      filename=FullModelTestingRunsDirName+"/run"+to_string(itest)+"/obs.txt";
      char obsnameread[200];;
      fptr_in=fopen(filename.c_str(),"r");
      if(fptr_in==NULL)
         printf("file not open?????\n");
      fscanf(fptr_in,"%s ",obsnameread);
      do{
         while(obsnameread[0]=='#'){
            fgets(cdummy,200,fptr_in);
            fscanf(fptr_in,"%s",obsnameread);
         }
         if(!feof(fptr_in)){
            fscanf(fptr_in,"%lf",&realY);
            if(!feof(fptr_in)){
               iY=observableinfo->GetIPosition(obsnameread);
               averagesigma2[iY]+=(realY-Y[iY])*(realY-Y[iY]);
               if(iY>=0 && iY<NObs)
                  fprintf(fptr_YvsY[iY],"%10.3e %10.3e %10.3e\n",realY,Y[iY],SigmaY_emulator[iY]);
               if(iY<NYY){
                  filename="smooth_data/outputstuff/fullmodel_testdata/YvsY_obs"+to_string(iY)+".txt";
                  fprintf(fptr_YvsY[iY],"%10.3e %10.3e %10.3e\n",realY,Y[iY],SigmaY_emulator[iY]);
               }
               if(fabs(realY-Y[iY])<SigmaY_emulator[iY])
                  nfit[iY]+=1;
            }
            //fgets(cdummy,200,fptr);
            fscanf(fptr_in,"%s",obsnameread);
         }
      }while(!feof(fptr_in));
      fclose(fptr_in);
      
   }
   for(iY=0;iY<NObs;iY++){
      fclose(fptr_YvsY[iY]);
   }
   
   for(iY=0;iY<NObs;iY++){
      averagesigma2[iY]=averagesigma2[iY]/double(ntest);
      CLog::Info("iY="+to_string(iY)+": <(Y-Yreal)^2>/SigmaY^2_emulator="+to_string(averagesigma2[iY]/(SigmaY_emulator[iY]*SigmaY_emulator[iY]))+"\n");
      fitpercentage=100.0*double(nfit[iY])/double(ntest);
      CLog::Info("percent < 1 sigma = "+to_string(fitpercentage)+"\n");
      ifit=floorl(fitpercentage);
      nfitpercent[ifit]+=1;
   }
   
   filename="smooth_data/output_stuff/fitpercentage.txt";
   FILE *fptr_percent=fopen(filename.c_str(),"w");
   for(ifit=0;ifit<100;ifit++){
      fprintf(fptr_percent,"%3d %g\n",ifit,double(nfitpercent[ifit]));
   }
   
   fclose(fptr_percent);

}
