#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>

using namespace std;

int main(int argc,char *argv[]){
	vector<int> list;
	int J;
	list.resize(5);
	list[0]=0; list[1]=1; list[2]=2; list[3]=3; list[4]=4;
	std::sort(list.begin(), list.end());
	J=0;
	do {
		J+=1;
		printf("%2d: ",J);
		for(int i=0;i<5;i++)
			printf("%2d ",list[i]);
		printf("\n");
	}while (std::next_permutation(list.begin(), list.end()));
}