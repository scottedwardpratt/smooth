#include "msu_smooth/simplex.h"
#include "msu_smooth/modelparinfo.h"
#include <cstdlib>

using namespace std;
using namespace NBandSmooth;
using namespace NMSUUtils;

CSimplexSampler::CSimplexSampler(){
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
	RGauss=parmap.getD("Simplex_RGauss",1.5);
	RGauss1=parmap.getD("Simplex_RGauss1",1.0);
	RGauss2=parmap.getD("Simplex_RGauss2",2.0);
}

void CSimplexSampler::SetThetaSimplex(){
	bool zeropoint=parmap.getB("Simplex_AddPointAtOrigin",false);
	if(TrainType==1){
		if(zeropoint)
			SetThetaType1Plus();
		else
			SetThetaType1();
	}
	else if(TrainType==2){
		if(zeropoint)
			SetThetaType2Plus();
		else
			SetThetaType2();
	}
	else if(TrainType==3){
		SetThetaType3();
	}
	else{
		CLog::Fatal("Inside CSimplexSampler::SetThetaSimplex, TrainType must be 1 or 2 or 3\n");
	}
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

	RTrain=fabs(RGauss);
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

void CSimplexSampler::SetThetaType1Plus(){
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

	RTrain=fabs(RGauss);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(priorinfo->type[ipar]=="gaussian")
				ThetaTrain[itrain][ipar]*=(RTrain/R);
			else{
				ThetaTrain[itrain][ipar]*=(RTrain/R);
			}
		}
	}
	
	// Add point at origin (differentiates this from type 1
	NTrainingPts+=1;
	ThetaTrain.resize(NTrainingPts);
	ThetaTrain[NTrainingPts-1].resize(NPars);
	for(ipar=0;ipar<NPars;ipar++)
		ThetaTrain[NTrainingPts-1][ipar]=0.0;
	
	
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
	unsigned int ipar,itrain,jtrain,N1,n,type2sign;
	double R,z,RTrain1,RTrain2;
	type2sign=parmap.getI("Simplex_Type2Sign",-1);
	if(NPars==2)
		type2sign=1;
	
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
	for(itrain=0;itrain<=NPars;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]=type2sign*ThetaTrain[itrain][ipar];
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
	RTrain1=fabs(RGauss1);
	RTrain2=fabs(RGauss2);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(itrain<=NPars){
				ThetaTrain[itrain][ipar]*=(RTrain1/R1);
			}
			else{
				ThetaTrain[itrain][ipar]*=(RTrain2/R2);
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

void CSimplexSampler::SetThetaType2Plus(){
	unsigned int ipar,itrain,jtrain,N1,n,type2sign=-1;
	double R,z,RTrain1,RTrain2;
	type2sign=parmap.getI("Simplex_Type2Sign",-1);
	
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
	
	for(itrain=0;itrain<=NPars;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]=type2sign*ThetaTrain[itrain][ipar];
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
	RTrain1=fabs(RGauss1);
	RTrain2=fabs(RGauss2);
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		for(ipar=0;ipar<NPars;ipar++){
			if(itrain<=NPars){
				ThetaTrain[itrain][ipar]*=(RTrain1/R1);
			}
			else{
				ThetaTrain[itrain][ipar]*=(RTrain2/R2);
			}
		}
	}
	
	// Add point at origin (differentiates this from type 2
	NTrainingPts+=1;
	ThetaTrain.resize(NTrainingPts);
	ThetaTrain[NTrainingPts-1].resize(NPars);
	for(ipar=0;ipar<NPars;ipar++)
		ThetaTrain[NTrainingPts-1][ipar]=0.0;
	

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

void CSimplexSampler::SetThetaType3(){
	unsigned int ipar,itrain,jtrain,N1,n;
	double R,z,RTrain;
	RTrain=fabs(RGauss);
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
	
	NTrainingPts=2*NPars+2;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=NTrainingPts/2;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=-ThetaTrain[itrain-NTrainingPts/2][ipar];
	}
	
	
	N1=NTrainingPts;
	n=N1;
	NTrainingPts+=(NPars+1)*NPars/2;
	ThetaTrain.resize(NTrainingPts);
	for(itrain=N1;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
	}
	for(itrain=1;itrain<NPars+1;itrain++){
		for(jtrain=0;jtrain<itrain;jtrain++){
			for(ipar=0;ipar<NPars;ipar++){
				ThetaTrain[n][ipar]=0.5*(ThetaTrain[itrain][ipar]+ThetaTrain[jtrain][ipar]);
			}
			n+=1;
		}
	}
	
	for(itrain=0;itrain<NTrainingPts;itrain++){
		R=0.0;
		for(ipar=0;ipar<NPars;ipar++){
			R+=ThetaTrain[itrain][ipar]*ThetaTrain[itrain][ipar];
		}
		R=sqrt(R);
		for(ipar=0;ipar<NPars;ipar++){
			ThetaTrain[itrain][ipar]*=(RTrain/R);
		}
	}
	
	// Add point at origin (differentiates this from type 2
	NTrainingPts+=1;
	ThetaTrain.resize(NTrainingPts);
	ThetaTrain[NTrainingPts-1].resize(NPars);
	for(ipar=0;ipar<NPars;ipar++)
		ThetaTrain[NTrainingPts-1][ipar]=0.0;
	

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

void CSimplexSampler::SetThetaTrain(vector<vector<double>> &theta){
	unsigned int itrain,ipar;
	for(itrain=0;itrain<ThetaTrain.size();itrain++)
		ThetaTrain[itrain].clear();
	ThetaTrain.clear();
	NTrainingPts=theta.size();
	ThetaTrain.resize(NTrainingPts);
	for(itrain=0;itrain<NTrainingPts;itrain++){
		ThetaTrain[itrain].resize(NPars);
		for(ipar=0;ipar<NPars;ipar++)
			ThetaTrain[itrain][ipar]=theta[itrain][ipar];
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

