// sphereDecoder.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <array>
#include <iostream>
#include "nrSphereDecorder.h"
#include "nrModulation.h"
#include "mwvc.h"
#include <random>
#include <numeric>
#include <algorithm>
#include <chrono>
#include </home/hliucp/Desktop/VertexCover/MySphereDecoder-CPP/eigen-3.4.0/Eigen/Dense>
#include </home/hliucp/Desktop/VertexCover/MySphereDecoder-CPP/eigen-3.4.0/Eigen/QR>
using namespace std;
//using namespace std::chrono; //::system_clock;

#include <sstream>
#include <fstream>

int binomialCoefficients(int n, int k);
Matrix2Dd example_run(int B, int n, int M, double D);
int main()    // int B, int M, int n  int argc, char *argv[]
{	
	ofstream file;
	//Veci Brange = {6,7};
	for (int M = 4; M<=6;M++){
	for (int B = 17; B < 31; B++) { // 2*M
		// int B=9; int n=9; int M=4; double D=sqrt(6-1e-6);    int B = Brange[i]; int M=4;
		// int B=7; int n=7; int M=3; int D = 4; double D1=sqrt(D-1e-6); cutoff_time = 1000;
		int n=B; int D = 4; double D1=sqrt(D-1e-6); cutoff_time = 1800;
		start = chrono::steady_clock::now();
		cout << "n=" << n << ", B=" << B << ", M=" << M << ", D=" << D;
		Matrix2Dd edgepairss = example_run(B, n, M, D1);
		double SD_comp_time = TimeElapsed();
		cout<<"SD_comp_time:" << SD_comp_time<< endl;
		bool a1= (TimeElapsed() >= cutoff_time);

		int remain1 = 0;
		uint seed = 5;
		for(int vcrun=0;vcrun<3;vcrun++){
			start = chrono::steady_clock::now();
			int v_numm = binomialCoefficients(B, M);
			BuildInstance(v_numm, edgepairss);
			if (mode < 0 || mode > 3)
			{
				mode = 0;
			}
			seed = seed+1;
			srand(seed);
			cutoff_time = 500;

			ConstructVC(); 
			LocalSearch();
			if (CheckSolution() == 1)
			{	
				int remain = v_num - best_weight; 
				cout << ", A =  " << remain << "; time: SD: " << SD_comp_time<< ", VC:" << best_comp_time << endl;
				if(remain1 < remain){
					string tablefilename = "/home/hliucp/Desktop/VertexCover/MySphereDecoder-CPP/table.txt";
					file.open(tablefilename, ios::app);  //open a file  VCforSD/
					file << endl;
					file << "n=" << n << ", B=" << B << ", M=" << M << ", D=" << D << endl;
					file << "total nodes: " << v_num << endl;
					file << "remained nodes: " << remain << endl;
					file << " time: SD: " << SD_comp_time << ", VC:" << best_comp_time << endl;
					file.close(); 

					Vecs pairlist = getNodeIdx(M, B);
					string filename = "/home/hliucp/Desktop/VertexCover/MySphereDecoder-CPP/B"+ to_string(B)+'M'+ to_string(M)+'D'+ to_string(D)+".txt";
					file.open(filename);
					file << endl;
					file << "remained nodes are ";
					file << remain << endl;
					for (int v = 1; v < v_num+1; v++)
					{	
						int a = *(best_v_in_c+v);
						if(a==0){
							// int sizepair = pairlist.size();
							string codewordd = pairlist[v-1];
							file<<codewordd<<",";
						}
					}
					file << endl;
					file.close();           //close it
					remain1 = remain;
				}

			}
			else
			{
				cout << ", the solution is wrong." << endl;
			}
			if (a1)
			{
				cout << "n=" << n << ", B=" << B << ", M=" << M << ", D=" << D<<"didn't finish in 1000s!"<<endl;
				break;
			}
		}
	}
	}
	FreeMemory();

	return 0;
}

int binomialCoefficients(int n, int k) {
   if (k == 0 || k == n)
   return 1;
   return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);
}

Matrix2Dd example_run(int nTxAnts = 5, int nRxAnts = 4, int M = 2, double D=2)
{
	// double N0 = 0.1; //noise
	// cout<<"N0 = "<< N0 <<endl;
	// vector<string> moduTypes{ "16qam","16qam","64qam" };
    // vector<string> moduTypes{ "qpsk","qpsk","qpsk","qpsk" };
	vector<string> moduTypes{ "3psk"};
	for(int i = 1;i<=nTxAnts-1;i++){
		moduTypes.push_back("3psk");
	}

	SphereDecoder sd(nRxAnts, nTxAnts, moduTypes);

	// random engines
	default_random_engine random_engine;
	bernoulli_distribution  bern_dist;
	normal_distribution<double> norm_dist(0, 1);

	int K = sd.getTotNumBits();
	Veci msg(nTxAnts); //msg.reserve(K);// generate msg
	// for (unsigned j = 0; j < K; j++)
	// 	msg.push_back(bern_dist(random_engine));
	// msg = {0,1,0,0,1,1,0,0,1,0};
	// msg = {1,1,1,1,1,1,1,1,1,1};
	Veci offsets = sd.getOffsets();
	Veci Ks = sd.getKs();

	ComplexVec txSymbs(nTxAnts); // txSymbs.reserve(nTxAnts);
	Veci bitlabels(nTxAnts);
	// cout<< "This is txSymbols:"<<endl;
	// for (auto i = 0; i < nTxAnts; i++) {
	// 	bitlabels = Veci(msg.begin() + offsets[i], msg.begin() + offsets[i] + Ks[i]);
	// 	// cout<< nrModuMapper(bitlabels, moduTypes[i])[0];
	// 	txSymbs.push_back(nrModuMapper(bitlabels, moduTypes[i],M)[0]);
	// }



	//flat channel
	random_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
	
	
	// using namespace Eigen;
	// Matrix<complex<double>, Dynamic, Dynamic> H;
	// H.setIdentity(nRxAnts, nTxAnts);
	ComplexMatrix2D H(nRxAnts);
	// Generate eye matrix
	for(int i = 0; i<=nRxAnts-1;i++){
		assert(nTxAnts==nRxAnts);
		ComplexVec Hrows(nTxAnts);
		Hrows[i]=1;
		H[i] = Hrows;
	}
	// int HH[4][5] = {{1,2,3,4,5},{5,4,3,2,1},{1,2,3,4,5},{5,4,3,2,1}};
	// for (auto i = 0; i < nRxAnts; i++) {
	// 	H[i].reserve(nTxAnts);
	// 	for (auto j = 0; j < nTxAnts; j++) {
	// 		H[i].push_back(1.0 / sqrt(2.0 * nTxAnts) * complex<double>(norm_dist(random_engine), norm_dist(random_engine)));	
	// 	}
	// }
	

	// rxSymbls
	ComplexVec rxSymbs; rxSymbs.reserve(nRxAnts);
	for (auto i = 0; i < nRxAnts; i++) {
		rxSymbs.push_back(0);
	}
	// decoded bits
	Matrix2Dd hardBits = sd.hardSphereDecode(H, rxSymbs,M,D), matrixVector = {};
	
	// Vecd dVector = {};
	// for(int i=0; i<=hardBits.size(); i++){
	// 	dVector.push_back(hardBits[i][0]);
	// }
	Veci V(hardBits.size());
	std::iota(V.begin(),V.end(),0); //Initializing

	std::sort( V.begin(),V.end(), [&](int i,int j){return hardBits[i]<hardBits[j];} );
	for(int i=0; i<hardBits.size(); i++){
		matrixVector.push_back(hardBits[V[i]]);
	}
	
	// return std::make_tuple(matrixVector,dVector);
	return matrixVector;
}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started:
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file