#include "msu_smooth/simplex.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CSimplexSampler::CSimplexSampler(){
	CparameterMap parmap;
	randy=new Crandy(123);
	parmap.ReadParsFromFile("smooth_data/smooth_parameters/simplex_parameters.txt");
	string logfilename=parmap.getS("Simplex_LogFileName","Screen");
	if(logfilename!="Screen"){
		CLog::Init(logfilename);
	}
	TrainType=parmap.getI("Simplex_TrainType",1);
	string prior_info_filename="smooth_data/Info/modelpar_info.txt";
	priorinfo=new CPriorInfo(prior_info_filename);
	CModelParameters::priorinfo=priorinfo;
	NPars=priorinfo->NModelPars;
	RGauss=parmap.getD("Simplex_RGauss",1.0);
	RGauss1=parmap.getD("Simplex_RGauss1",0.5);
	RGauss2=parmap.getD("Simplex_RGauss2",1.0);
}

void CSimplexSampler::SetThetaSimplex(){
	if(TrainType==1)
		SetThetaType1();
	else if(TrainType==2)
		SetThetaType2();
	else{
		CLog::Fatal("Inside CSimplexSampler::SetThetaSimplex, TrainType must be 1 or 2\n");
	}
	CLog::Info("NTrainingPts="+to_string(NTrainingPts)+"\n");
}

void CSimplexSampler::SetThetaType1(){
	unsigned int ipar,itrain,jtrain;
	double R,z,RTrain;
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}

	RTrain=RGauss*sqrt(double(NPars)/3.0);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="gaussian")
				ThetaTrain[itrain][ipar]*=(RTrain/R);
			else{
				ThetaTrain[itrain][ipar]*=(RTrain/R);
			}
		}
	}

	double Rmax=0.95;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="uniform"){
				if(fabs(ThetaTrain[itrain][ipar])>Rmax){
					ThetaTrain[itrain][ipar]*=(Rmax/fabs(ThetaTrain[itrain][ipar]));
				}
			}
		}
	}
	
	
}

void CSimplexSampler::SetThetaType2(){
	unsigned int ipar,itrain,jtrain,N1,n;
	double R,z,RTrain1,RTrain2;
	RTrain2=sqrt(NPars);
	RTrain1=RTrain2/2.0;
	
	NTrainingPts=NPars+1;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=0.0;
	}
	R=1.0;
	ThetaTrain[0][0]=-R;
	ThetaTrain[1][0]=R;
	for(itrain=2;itrain<NTrainingPts;itrain++){
		z=R*itrain/sqrt(double(itrain*itrain)-1.0);
		for(jtrain=0;jtrain<itrain;jtrain++){
			ThetaTrain[jtrain][itrain-1]=-z/double(itrain);
		}
		ThetaTrain[itrain][itrain-1]=z;
		R=z;
	}
	
	N1=NTrainingPts;
	n=N1;
	NTrainingPts+=N1*(N1-1)/2;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=N1;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
	}
	for(itrain=1;itrain<N1;itrain++){
		for(jtrain=0;jtrain<itrain;jtrain++){
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[n][ipar]=0.5*(ThetaTrain[itrain][ipar]+ThetaTrain[jtrain][ipar]);
			}
			n+=1;
		}
	}
	
	double R1,R2;
	itrain=0;
	R1=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		R1+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
	R1=sqrt(R1);
	itrain=NTrainingPts-1;
	R2=0.0;
	for(ipar=0;ipar<NPars;ipar++)
		R2+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
	R2=sqrt(R2);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(itrain<=NPars){
				ThetaTrain[itrain][ipar]*=(RGauss1/R1);
			}
			else{
				ThetaTrain[itrain][ipar]*=(RGauss2/R2);
			}
		}
	}


	double Rmax=0.95;
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="uniform"){
				if(fabs(ThetaTrain[itrain][ipar])>Rmax){
					ThetaTrain[itrain][ipar]*=(Rmax/fabs(ThetaTrain[itrain][ipar]));
				}
			}
		}
	}
	
}

void CSimplexSampler::WriteModelPars(){
	FILE *fptr;
	string filename,dirname,command;
	unsigned int itrain,ipar;
	vector<CModelParameters *> modelparameters(NTrainingPts);

	for(itrain=0;itrain<NTrainingPts;itrain++){
		modelparameters[itrain]=new CModelParameters();
		for(ipar=0;ipar<NPars;ipar++){
			modelparameters[itrain]->Theta[ipar]=ThetaTrain[itrain][ipar];
		}
		modelparameters[itrain]->TranslateTheta_to_X();
	}
	for(itrain=0;itrain<NTrainingPts;itrain++){
		dirname="smooth_data/modelruns/run"+to_string(itrain);
		command="mkdir -p "+dirname;
		system(command.c_str());
		filename=dirname+"/mod_parameters.txt";
		fptr=fopen(filename.c_str(),"w");
		for(ipar=0;ipar<NPars;ipar++){
			fprintf(fptr,"%s %g\n",
			priorinfo->parname[ipar].c_str(),modelparameters[itrain]->X[ipar]);
		}
		fclose(fptr);
	}
	
}
