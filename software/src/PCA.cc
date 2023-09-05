#include "msu_smooth/PCA.h"

using namespace std;

PCA::PCA(string filename){
	CparameterMap *parmap=new CparameterMap();
	parmap->ReadParsFromFile(filename);

	observable_info=new CObservableInfo(parmap->getS("SmoothEmulator_ObservableInfoFileName","Info/observable_info.txt"));
	


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
	
	A = new Eigen::MatrixXd(Nobs,Nobs);
	A->setZero();
	
	for(iy=0;iy<Nobs;iy++){
		for(jy=0;jy<Nobs;jy++){
			for(irun=0;irun<nruns;irun++){
				(*A)(iy,jy)+=(Ytilde[irun][iy]*Ytilde[irun][jy])/double(nruns);
			}
		}
	}
				
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(*A);
	eigvals = es.eigenvalues();
	eigvecs = es.eigenvectors();

	ofstream PCAInfoFile("PCA_Info/info.txt");
	
	if (PCAInfoFile.is_open()){
		PCAInfoFile << "Nobs  Nruns" << endl;
		PCAInfoFile << Nobs <<  " " << nruns << endl;
		PCAInfoFile << "<Y>" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << Ybar[i] " ";
		PCAInfo/file << endl;
		PCAInfoFile << "SigmaY" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << SigmaY[i] " ";
		PCAInfo/file << endl;
		PCAInfoFile << "EigenValues" << endl;
		for(iy=0;iy<Nobs;iy++)
			PCAInfoFile << eigvals(iy) " ";
		PCAInfo/file << endl;
		PCAInfoFile << "EigenVectors" << endl;
		for(iy=0;iy<Nobs;iy++){
			for(jy=0;jy<Nobs;jy++){
				PCAInfoFile << eigvecs(iy,jy) " ";
			}
			PCAInfo/file << endl;
		}
	}
	PCAInfoFile.close();

}

void PCA::ReadPCA(){
	ifstream eigvalFile("PCA_Info/eigvals.txt");
	string line;
	vector<double> v;

	while (std::getline(eigvalFile, line)){
		v.push_back(stod(line));
	}


	int matrix_size = v.size();

	eigvals.resize(matrix_size, matrix_size);
	eigvals.setZero();



	for (size_t i = matrix_size; i > 0; i--) {
		eigvals(i-1,i-1) = v[matrix_size - i];
	}

	eigvalFile.close();

	ifstream eigvecFile("PCA_Info/eigvecs.txt");
	eigvecs.resize(matrix_size, matrix_size);
	int col = 0;

	while (std::getline(eigvecFile, line)) {
		int row = 0;
		stringstream str_strm;
		str_strm << line;
		string temp_str;
		double temp_dob;
		while(!str_strm.eof()){
			str_strm >> temp_str;
			if (stringstream(temp_str) >> temp_dob) {
				eigvecs(col , matrix_size - 1 - row) = temp_dob;
				row++;
			}
			temp_str = "";
		}
		col++;
	}
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
