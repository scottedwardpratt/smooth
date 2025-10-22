#ifndef __SMOOTH_MASTER_H__
#define __SMOOTH_MASTER_H__
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <sstream>
#include "msu_smoothutils/parametermap.h"
#include "msu_smoothutils/misc.h"
#include "msu_smoothutils/randy.h"
#include "msu_smooth/pca.h"
#include "msu_smooth/emulator.h"
#include "msu_smooth/modelparinfo.h"
#include "msu_smooth/smooth.h"
#include "msu_smoothutils/log.h"
#include "msu_smooth/observableinfo.h"
#include "msu_smooth/priorinfo.h"
#include "msu_smooth/traininginfo.h"

using namespace NMSUUtils;

namespace NBandSmooth{
   class CSmoothEmulator;
   class CTrainingInfo;
   class CPriorInfo;
   class CObservableInfo;
   class CModelParInfo;
   
   class CSmoothMaster{
   public:
      bool UsePCA;
      unsigned int TrainType;
      CSmoothMaster();
      CparameterMap *parmap;
      unsigned int NPars;
      vector<CSmoothEmulator *> emulator;
      CTrainingInfo *traininginfo;
      CObservableInfo *observableinfo;
      CPriorInfo *priorinfo;
      string FullModelRunDirName;
      Crandy *randy;
      CSmooth *smooth;
      string CoefficientsDirName,SurmiseTrainingPointsFileName,SurmixeTrainingObsFileName;
      string SmoothEmulator_TrainingFormat;
      double pca_minvariance,fitpercentage;
      vector<bool> pca_ignore;
      int GetNPars(){
         return NPars;
      }
      int GetNObs(){
         return observableinfo->NObservables;
      }
      
      void ReadTrainingInfo();
      //void GenerateCoefficientSamples();
      void CalcAllSigmaALambda();
      void TuneAllY(); // tune all observables
      void TuneY(string obsname); // tune one observable
      void TuneY(unsigned int iY); // tune one observable
      void TuneAllY(double LAMBDA); // tune all observables
      void TuneY(string obsname,double LAMBDA); // tune one observable
      void TuneY(unsigned int iY,double LAMBDA); // tune one observable
      
      void GetAllY(CModelParameters *modelpars,vector<double> &Y,vector<double> &SigmaY_emulator);
      void GetAllY(vector<double> &theta,vector<double> &Y,vector<double> &SigmaY_emulator);
      void GetAllYOnly(CModelParameters *modelpars,vector<double> &Y);
      void GetAllYOnly(vector<double> &theta,vector<double> &Y);
      
      
      void GetY(unsigned int iY,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
      void GetY(unsigned int iY,vector<double> &theta,double &Y,double &SigmaY_emulator);
      void GetY(string obsname,CModelParameters *modelpars,double &Y,double &SigmaY_emulator);
      void GetY(string obsname,vector<double> &theta,double &Y,double &SigmaY_emulator);
      
      double GetYOnly(unsigned int iY,CModelParameters *modelpars);
      double GetYOnly(unsigned int iY,vector<double> &theta);
      double GetYOnly(string obsname,CModelParameters *modelpars);
      double GetYOnly(string obsname,vector<double> &theta);
      double GetYOnly(int iY,vector<double> theta);
      double GetYOnlyPython(int DiY,vector<double> theta);
      vector<double> GetYSigmaPython(int DiY,vector<double> theta);
      
      double GetUncertainty(string obsname,vector<double> &Theta);
      double GetUncertainty(unsigned int iY,vector<double> &theta);
      double GetUncertainty(int iY,vector<double> theta);
      
      void TestAtTrainingPts();
      void TestAtTrainingPts(string obsname);
      void TestAtTrainingPts(unsigned int iY);
      void TestVsFullModelAlt();
      void TestVsFullModel();
      
      vector<double> GetXFromTheta(vector<double> Theta);
      vector<double> GetThetaFromX(vector<double> X);
      
   };
   
};

//#include "msu_smooth/smoothbind.h"

#endif
