#include "smooth.h"
using namespace std;

CEmulator::CEmulator(unsigned int NPars_set,unsigned int MaxRank_set,unsigned int NSmooth_set){
	NSmooth=NSmooth_set;
	MaxRank=MaxRank_set;
	NPars=Npars_set;
	NTrainingPts=0;
	smooth.resize(NSmooth);
	for(unsigned int ismooth=0;ismooth<NSmooth;ismooth++){
		smooth[ismooth]=new CSmooth(this);
	}
}

void SetNTrainingPts(unsigned int NTrainingPts_set);{
	if(NTrainingPts!=NTrainingPts_set){
		delete gslmatrix;
	}
	NTrainingPts=NTrainingPts_set;
	gslmatrix=new CGSLMatrix_Real(NTrainingPts);
}

void CEmulator::SolveForSmoothPars(){
	//solve for AAx=y, where AA=theta^m... , AA are first NTrainingPts values of theta^m and x are The A coefficients and y are the training values
	
	vector<double> xvalues;(NTrainingPts);
	vector<double> yvalues;(NTrainingPts);
	vector<vector<double>> AA;
	xvalues.resize(NTrainingPoints);
	yvalues.resize(NTrainingPoints);
	AA.resize(NTrainingPts);
	for(unsigned int it=0;it<NTrainingPts;it++){
		AA[it].resize(NTrainingPts);
		for(unsigned in jt=0;jt<NTrainingPts;jt++)
			A[it][jt]=0.0;
	}
	for(unsigned int itheta=0;itheta<NTrainingPts;itheta++){
		smooth.GetPrefactors(NTrainingPts,trainingthetas[itheta],AA[itheta]);
	}
	gslmatrix->SolveLinearEqs(trainingvalues,AA,x);
	smooth.SetParsFromX(x);
	
}
