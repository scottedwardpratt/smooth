#include "msu_smooth/PCA.h"

using namespace std;

PCA::PCA(string filename){
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(filename);

	observable_info=new CObservableInfo("Info/observable_info.txt");

	modelruns_dirname=parmap->getS("SmoothEmulator_ModelRunDirName","modelruns");
	string NTrainingStr = parmap->getS("SmoothEmulator_TrainingPts","1");
	
	stringstream ss(NTrainingStr);
	string token;

	while(getline(ss, token, ',')) {
		size_t pos = token.find("-");
		if (pos != string::npos) {

			int start = stoi(token.substr(0, pos));
			int end = stoi(token.substr(pos+1));

			for (int i = start; i <= end; i++)
				NTrainingList.push_back(i);
		}
		else {

			NTrainingList.push_back(stoi(token));
		}
	}
	nruns=NTrainingList.size();

	Y.resize(nruns);
	SigmaY.resize(nruns);

}

void PCA::CalcPCA(){
	vector<vector<double>> Ytilde,Y;
	vector<double> Ytildebar;
	char filename[300],obs_name[300];
	double y,sigmay;
	int iy,jy,irun,nruns=NTrainingList.size();
	Nobs=observable_info->NObservables;
	Eigen::MatrixXd *A;
	FILE *fptr;
	Y.resize(nruns);
	Ytilde.resize(nruns);
	for(irun=0;irun<nruns;irun++){
		Y[irun].resize(Nobs);
		Ytilde[irun].resize(Nobs);
	}
	Ybar.resize(Nobs);
	SigmaY.resize(Nobs);
	for(iy=0;iy<Nobs;iy++)
		SigmaY[iy]=Ybar[iy]=0.0;
	
	
	for(irun=0;irun<nruns;irun++){
		snprintf(filename,300,"%s/run%d/obs.txt",modelruns_dirname.c_str(),NTrainingList[irun]);
		fptr=fopen(filename,"r");
		while(!feof(fptr)){
			fscanf(fptr,"%s %lf %lf",obs_name, &y, &sigmay);
			if(!feof(fptr)){
				iy=observable_info->GetIPosition(obs_name);
				SigmaY[iy]+=sigmay/double(nruns);
				Y[irun][iy]=y;
				Ybar[iy]+=y/double(nruns);
			}
		}
		fclose(fptr);
	}
	
	for(iy=0;iy<Nobs;iy++){
		for(irun=0;irun<nruns;irun++){
			Ytilde[irun][iy]=(Ytilde[irun][iy]-Ybar[iy])/SigmaY[iy];
		}
	}
	
	// Calculate Covariance Matrix
	
	A = new Eigen::MatrixXd(Nobs,Nobs);
	A->setZero();
	
	for(iy=0;iy<Nobs;iy++){
		for(jy=0;jy<Nobs;jy++){
			for(irun=0;irun<nruns;irun++){
				(*A)(iy,jy)+=(Ytilde[irun][iy]*Ytilde[irun][jy])/double(nruns);
			}
		}
	}
	
	// Write Transformation Info
				
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*A);
	eigvals = es.eigenvalues();
	eigvecs = es.eigenvectors();

	ofstream PCAInfoFile("PCA_Info/transformation_info.txt");
	
	if (PCAInfoFile.is_open()){
		PCAInfoFile << "Nobs  Nruns" << endl;
		PCAInfoFile << Nobs <<  " " << nruns << endl;
		PCAInfoFile << "<Y>" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << Ybar[iy] << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "SigmaY" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << SigmaY[iy] << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "EigenValues" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << eigvals(iy) << " ";
		PCAInfoFile << endl;
		PCAInfoFile << "EigenVectors" << endl;
		for(iy=0;iy<Nobs;iy++){
			for(jy=0;jy<Nobs;jy++){
				PCAInfoFile << eigvecs(iy,jy) << " ";
			}
			PCAInfoFile << endl;
		}
	}
	PCAInfoFile.close();
	
	delete A;
	
	// Write PCA Observables for Training Pts
	
	vector<vector<double>> Z;
	int iz;
	string pcaname;
	
	Z.resize(nruns);
	for(irun=0;irun<nruns;irun++){
		Z[irun].resize(nruns);
	}
	
	for(irun=0;irun<nruns;irun++){
		for(iz=0;iz<Nobs;iz++){
			Z[irun][iz]=0.0;
			for(iy=0;iy<Nobs;iy++){
				Z[irun][iz]+=eigvecs(iz,iy)*Ytilde[irun][iy];
			}
		}
		snprintf(filename,300,"%s/run%d/pca_obs.txt",modelruns_dirname.c_str(),NTrainingList[irun]);
		fptr=fopen(filename,"w");
		for(iz=0;iz<Nobs;iz++){
			pcaname="z"+to_string(iz);
			fprintf(fptr,"%s  %g  0.0z\n",pcaname.c_str(),Z[irun][iz]);	
		}
		fclose(fptr);
	}
	
	// Write Info File for PCA
	
	double SA0Zsquared,SA0Y,SA0Z;
	
	fptr=fopen("Info/pca_info.txt","w");
	
	for(int iz=0;iz<Nobs;iz++){
		SA0Zsquared=0.0;
		for(iy=0;iy<nruns;iy++){
			SA0Y=observable_info->SigmaA0[iy];
			SA0Zsquared+=SA0Y*SA0Y*eigvecs(iz,iy)*eigvecs(iz,iy);
		}
		SA0Z=sqrt(SA0Zsquared);
		pcaname="z"+to_string(iz);
		fprintf(fptr,"%s %g\n",pcaname.c_str(),SA0Z);
	}
	
	fclose(fptr);
	

}

void PCA::ReadPCATransformationInfo(){
	char dummy[200];
	int iy,jy;
	FILE *fptr=fopen("PCA_Info/transformation_info.txt","r");
	fgets(dummy,200,fptr);
	fscanf(fptr,"%d %d",&Nobs,&nruns);
	Ybar.resize(Nobs);
	SigmaY.resize(Nobs);
	eigvals.resize(Nobs);
	eigvecs.resize(Nobs,Nobs);
	
	for(iy=0;iy<Nobs;iy++){
		fscanf(fptr,"%lf ",&Ybar[iy]);
	}
	fgets(dummy,200,fptr);
	fgets(dummy,200,fptr);
	for(iy=0;iy<Nobs;iy++){
		fscanf(fptr,"%lf ",&SigmaY[iy]);
	}
	fgets(dummy,200,fptr);
	fgets(dummy,200,fptr);
	for(iy=0;iy<Nobs;iy++){
		fscanf(fptr,"%lf",&eigvals(iy));
	}
	fgets(dummy,200,fptr);
	fgets(dummy,200,fptr);
	for(iy=0;iy<Nobs;iy++){
		for(jy=0;jy<Nobs;jy++){
			fscanf(fptr,"%lf ",&eigvecs(iy,jy));
		}
	}
	fclose(fptr);
	
}

void PCA::WriteZTraining(){
	int ifile;

	string command="rm -r -f PCA_Info/run*";
	system(command.c_str());

	for (int irun = 0; irun < nruns; irun++) {
		string filename;
		FILE *fptr;

		string dirname= "PCA_Info/run" +to_string(irun);
		command= "mkdir " +dirname;
		system(command.c_str());

		Eigen::VectorXd obsVec(Y[0].size());

		for (size_t i = 0; i < Y[0].size(); i++){
			obsVec(i) = Y[irun][i];
		}

		Eigen::VectorXd PCAVec(Y[0].size());

		PCAVec = eigvecs * obsVec;
		
		ifile=NTrainingList[irun];

		filename=modelruns_dirname+"/run"+to_string(ifile)+"/obs_pca.txt";
		fptr=fopen(filename.c_str(),"w");

		for(size_t i = 0; i < Y[0].size(); i++){
			string obsname = "pca" + to_string(i);
			fprintf(fptr,"%s %s %g\n","dummy" ,obsname.c_str(),PCAVec(i));
		}
		fclose(fptr);
	}
}

void PCA::TranslateZtoY(vector<double> &Z,vector<double> &Y,vector<double> &SigmaZ_emulator,vector<vector<double>> &SigmaY_emulator){
	int iy,jy,ky;
	SigmaY_emulator.resize(Nobs);
	for(iy=0;iy<Nobs;iy++){
		SigmaY_emulator.resize(Nobs);
	}
	for(iy=0;iy<Nobs;iy++){
		Y[iy]=0.0;
		for(ky=0;ky<Nobs;ky++){
			Y[iy]+=Z[ky]*eigvecs(ky,iy);
		}
		for(jy=0;jy<Nobs;jy++){
			SigmaY_emulator[iy][jy]=0.0;
			for(ky=0;ky<Nobs;ky++){
				SigmaY_emulator[iy][jy]+=SigmaZ_emulator[ky]*SigmaZ_emulator[ky]*eigvecs(ky,iy)*eigvecs(ky,jy);
			}
			SigmaY_emulator[iy][jy]*=SigmaY[iy]*SigmaY[jy];
		}
	}
	
}
