//--------------------------------------------------------------------------------------------------------------
//                                       ** sphere decoder **
//						Author : Dr J Mao  Email : juquan.justin.mao@gmail.com  2021 Dec
//          source paper :  "Soft-Output Sphere Decoding: Performance and Implementation Aspects"
//-------------------------------------------------------------------------------------------------------------
#include <iostream>
#include "nrSphereDecorder.h"
// #include "mwvc.h"
#include <numeric>
#include <algorithm>
#include <array>
#include </home/hliucp/Desktop/VertexCover/MySphereDecoder-CPP/eigen-3.4.0/Eigen/Dense>
#include </home/hliucp/Desktop/VertexCover/MySphereDecoder-CPP/eigen-3.4.0/Eigen/QR>
#include "assert.h" 

using namespace std;

SphereDecoder::SphereDecoder(int nRxAnts, int nTxAnts, vector<string>& moduTypes, double searchRadius)
{
	NumRxAnts = nRxAnts; // number of rx antennas
	NumTxAnts = nTxAnts; // number of transmit antennas

	// Search Radius controls tree pruning for counter hypothesis, 0: hard decision,inf: max-log-map
	// values in between is the tradeoff between complexityand performance, 0.2 is a thrumb of rule value.
	SearchRadius = searchRadius;

	//modulation handling
	if (moduTypes.size() != NumTxAnts) {
		if (moduTypes.size() == 1) {
			//string tmp = moduTypes[0];
			moduTypes.insert(moduTypes.end(), NumTxAnts - 1, moduTypes[0]);
		}
		else {
			std::cerr << "The number of modulations must be 1 or equal to number of Tx antennas!" << std::endl;
			abort();
		}
	}

	// constellations for each Tx Ants
	Cstltns.reserve(NumTxAnts);
	for (auto i = 0; i < NumTxAnts; i++) Cstltns.push_back(Constellation{ moduTypes[i] });

	// record number of bits per decoded symbol
	Ks.reserve(nTxAnts);
	for (auto i = 0; i < nTxAnts; i++) Ks.push_back(Cstltns[i].K);

	// total number of bits corresponding the decoded symbols
	Offsets.reserve(nTxAnts); // record the starting position for bit labels of each decoded symbols
	NumTotBits = 0;
	for (auto i = 0; i < nTxAnts; i++) {
		Offsets.push_back(NumTotBits);
		NumTotBits += Ks[i];
	}
}
SphereDecoder::~SphereDecoder()
{
}

void SphereDecoder::nodeInit(TreeNode& node)
{
	// initilize to root node;
	// add the first symbol of the root costellation
	node.level = NumTxAnts - 1;
	Constellation& cstltn = Cstltns[node.level];
	node.psv.push_front(cstltn.ConstlSymbs[0]);
	node.symIndices.push_front(0);
	node.bitLabels.insert(node.bitLabels.end(), cstltn.BitLabels[0].begin(), cstltn.BitLabels[0].end());
}

Vecd SphereDecoder::operator()(ComplexMatrix2D& H, ComplexVec& rxSymbs)
{
	//-----------------------------------------------------------------------------------------------
	// softbits = nrSphereDecoder(H, rxSymbs) uses sphere decoding algorithms (single tree seach)
	// for seeking the maximum - likelihood solution for a set of symbols transmitted over the MIMO channel.
	// ***Input***
	//	* H - complex channel matrix Nr-by-Nt, Nr is the number of rx antennas, Nt is the number of Tx antennas
	//	* rxSymbs - complex received symbol vector  Nr-by-1
	// ***Output*** softbits(llr) - soft bits (not sacled by 1/N0)
	//--------------------------------------------------------------------------------------------------

	// hypothesis maximum likelihood (ML) bits
	Veci bitLabelsML; bitLabelsML.reserve(NumTotBits);
	for (auto i = 0; i < NumTxAnts; i++)
		bitLabelsML.insert(bitLabelsML.end(), Cstltns[i].BitLabels[0].begin(), Cstltns[i].BitLabels[0].end());

	// QR Decomposition
	ComplexMatrix2D Qh, R; qrDecomposeEigen(H, Qh, R);

	// perform Q'*y, y is rx symbol vector
	ComplexVec Qhy; Qhy.reserve(NumTxAnts);
	for (auto i = 0; i < NumTxAnts; i++) {
		Qhy.push_back(inner_product(Qh[i].begin(), Qh[i].end(), rxSymbs.begin(), complex<double>(0.0)));
	}

	// initalize node to root
	TreeNode node; nodeInit(node);

	double di; //distance increment
	Vecd ped(NumTxAnts, 0.0);// partial euclidean distance
	//hypothesis minimum distances for ML bitsand its counterpart
	double lambdaML = std::numeric_limits<double>::infinity(); // minimum distance
	Vecd lambdaMLbar(NumTotBits, std::numeric_limits<double>::infinity());

	// varibles related to tree traveral
	bool isTravelDone = false;
	NumNodesVisted = 0; // only count leaf nodes

	while (!isTravelDone)
	{
		//----------------- disply node ----------------------------------------------------
		//cout << "node ["; for (auto e : node.symIndices) std::cout << e << " "; cout << "\b]" << endl;
		//
		// compute distance increment of current node
		di = norm(Qhy[node.level] -
			inner_product(node.psv.begin(), node.psv.end(), R[node.level].begin() + node.level, complex<double>(0.0)));

		// compute partial Euclidean distance, none root node distance increamented di from its parent node.
		ped[node.level] = (node.level == NumTxAnts - 1) ? di : ped[node.level + 1] + di;

		if (node.level == 0) { // leaf node
			NumNodesVisted += 1;
			if (ped.front() < lambdaML) { // smaller euclidean distance found
				//update the counter hypotheses
				for (auto i = 0; i < NumTotBits; i++) {
					if (node.bitLabels[i] != bitLabelsML[i]) lambdaMLbar[i] = lambdaML;
				}

				//update the hypotheses
				lambdaML = ped.front(); bitLabelsML = Veci(node.bitLabels.begin(), node.bitLabels.end());

				//LLR clipping for counter hypothese
				for (auto& e : lambdaMLbar) e = min(e, SearchRadius + lambdaML);
			}
			else { // No smaller euclidean distance found in the leaf node
				// update the counter hypotheses only
				for (auto i = 0; i < NumTotBits; i++) {
					if ((node.bitLabels[i] != bitLabelsML[i]) & (lambdaMLbar[i] > ped.front())) lambdaMLbar[i] = ped.front();
				}
			}
			isTravelDone = moveRight(node);
		}
		else { // none leaf node
			if (ped[node.level] < lambdaML)
				moveDown(node);
			else { // for hypothesis there is no need to go down
				//Check if there is a smaller Euclidean distance for the couther hypothesis in the sub - tree

				//  find the maximum lambadaMLBar for those bits not yet traversed
				int nBitsNotVisted = accumulate(Ks.begin(), Ks.begin() + node.level, 0);
				double lambaMax = *max_element(lambdaMLbar.begin(), lambdaMLbar.begin() + nBitsNotVisted);

				// find the maximum lambadaMLBar for those bits already traversed
				for (auto i = 0; i < NumTotBits - nBitsNotVisted; i++) {
					if (node.bitLabels[i] != bitLabelsML[nBitsNotVisted + i])
						if (lambdaMLbar[nBitsNotVisted + i] > lambaMax)
							lambaMax = lambdaMLbar[nBitsNotVisted + i];
				}

				if (ped[node.level] < lambaMax) // not punning as counter hypothesis  has a chance to update
					moveDown(node);
				else
					isTravelDone = moveRight(node); // prunning
			}
		}
	}
	// comput soft bits
	Vecd llr; llr.reserve(NumTotBits);
	for (auto i = 0; i < NumTotBits; i++) {
		llr.push_back((bitLabelsML[i] == 1) ? lambdaML - lambdaMLbar[i] : lambdaMLbar[i] - lambdaML);
	}
	return llr;
}

void combinate(int iPos, int iProc, int iTol, int iKey, Veci data,Veci des,Matrix2Di& per)
{   
	if(iProc > iTol)
	{//cout<<"51"<<endl;
		return;
	}
	if(iPos == iKey)
	{ // cout<<"51"<<endl;
	    per.push_back(des);
		// for(int i = 0;i < iKey; i++)
		// {//cout<<"31"<<endl;
		// 	cout<<des[i]<<" ";
		// }
		// cout<<endl;
		return;
	}
	else
	{  // cout<<"1"<<endl;
	// cout<<"31"<<endl;
		combinate(iPos,iProc+1,iTol,iKey,data,des,per);
		//cout<<"41"<<endl;
		des[iPos] = data[iProc];
		combinate(iPos+1,iProc+1,iTol,iKey,data,des,per);
	}
}

void nchoosek(Veci data, int M, Matrix2Di& per){
    // nchoosek for vector data
	int B = data.size();
    Veci temp(M);
    //cout<<"21"<<endl;
    combinate(0, 0, B, M, data, temp, per);
    //cout<<"1"<<endl;
    return;
} 

Veci findPrime(int primesize) {
	const int MAXN = 10000; 
	Veci prime;
	bool p[MAXN] = {0};	// 初始化，默认所有的数都是素数
	for (int i = 2; i < MAXN; ++i) {
		if (p[i] == false) {
			prime.emplace_back(i);
			if (prime.size() >= primesize) break;
			for (int j = i + i; j < MAXN; j += i) {
				// 将 i 的所有倍数标记为true，表示其不是素数
				p[j] = true;
			}
		}
	}
	return prime;
}

Vecs getNodeIdx(int M, int B){
	Matrix2Di idx1;
	Vecs pairlist;
	Veci pair1;
	Veci data(B);
	iota(data.begin(),data.end(),0);
	
	nchoosek(data, M, idx1);
	for(Veci ii:idx1){
		std::string idx1_s = {};
		pair1.assign(ii.begin(),ii.end()); 
		std::sort(pair1.begin(),pair1.end()); 
		for(int s:pair1) idx1_s.append(to_string(s)).append(1,'_');
		pairlist.push_back(idx1_s);
		}
	return pairlist; 
}

Matrix2Ds SphereDecoder::getNodesPairIdx(Veci E, int M){
	Veci idx1, idx_1, idx0;
	int m1=0, m_1=0;
	for (int i=0;i<E.size();i++){
		switch(E[i]){
			case 1: idx1.push_back(i); m1 += 1; break;
			case 0: idx0.push_back(i);break;
			case -1: idx_1.push_back(i); m_1 +=1;break;
		}
	}
	assert(m1 == m_1);
	Matrix2Di idxPair1;
	Matrix2Ds pairIdx,primePair1;
	Veci prime = findPrime(E.size()),pair1,pair_1;
	Vecd prime1,prime_1;
	std::string idx1_s = {},idx_1_s = {};
	// for(int i:idx1) prime1.push_back({1.0/prime[i]});
	// for(int i:idx_1) prime_1.push_back({1.0/prime[i]});
	// double aa1 =  accumulate(prime1.begin(), prime1.end(), 0);
	// double aa_1 = accumulate(prime_1.begin(), prime_1.end(), 0);
	switch (M-m1){
		case 0:
		for(int s:idx1) idx1_s.append(to_string(s)).append(1,'_');
		for(int s:idx_1) idx_1_s.append(to_string(s)).append(1,'_');
		pairIdx = {{idx1_s,idx_1_s}}; // idxPair1 = {{0}};
		break;
		case 1:   // right???
		for(int i:idx0){
			idx1_s = {},idx_1_s = {};
			pair1.assign(idx1.begin(),idx1.end()); pair_1.assign(idx_1.begin(),idx_1.end());
			pair1.push_back(i); std::sort(pair1.begin(),pair1.end());
			pair_1.push_back(i);std::sort(pair_1.begin(),pair_1.end());
			for(int s:pair1) idx1_s.append(to_string(s)).append(1,'_');
			for(int s:pair_1) idx_1_s.append(to_string(s)).append(1,'_');
			pairIdx.push_back({idx1_s,idx_1_s});
		} 
		break;
		default : 
		nchoosek(idx0, M-m1, idxPair1);   // right???
		for(Veci ii:idxPair1){
			pair1.assign(idx1.begin(),idx1.end()); pair_1.assign(idx_1.begin(),idx_1.end());
			for(int i:ii){
				pair1.push_back(i); 
				pair_1.push_back(i);
			}
			std::sort(pair1.begin(),pair1.end()); std::sort(pair_1.begin(),pair_1.end());
			idx1_s = {},idx_1_s = {};
			for(int s:pair1) idx1_s.append(to_string(s)).append(1,'_');
			for(int s:pair_1) idx_1_s.append(to_string(s)).append(1,'_');
			pairIdx.push_back({idx1_s,idx_1_s});
		}
	// 	for(int i=0;i<idxPair1.size();i++){
	// 		for(int j=0;j<idxPair1[i].size();j++){
	// 			int k = idxPair1[i][j];
	// 			primePair1[i].push_back({1.0/prime[k]});
	// 		}
	// 		double a = accumulate(primePair1[i].begin(), primePair1[i].end(), 0);
	// 		pairIdx.push_back({a + aa1,a + aa_1});
	// 	} 
	}
	return pairIdx;

}


Matrix2Dd SphereDecoder::hardSphereDecode(ComplexMatrix2D& H, ComplexVec& rxSymbs, int M, double SearchRadius)
{	
	//-----------------------------------------------------------------------------------------------
	// ***Input***
	//	* H - complex channel matrix Nr-by-Nt, Nr is the number of rx antennas, Nt is the number of Tx antennas
	//	* rxSymbs - complex received symbol vector  Nr-by-1
	// ***Output***
	//  * hard detected bits
	//--------------------------------------------------------------------------------------------------

	// hypothesis maximum likelihood (ML) bits
	Veci bitLabelsML; bitLabelsML.reserve(NumTotBits);
	for (auto i = 0; i < NumTxAnts; i++)
		bitLabelsML.insert(bitLabelsML.end(), Cstltns[i].BitLabels[0].begin(), Cstltns[i].BitLabels[0].end());

	ComplexMatrix2D Qh, R, Hh;
	ComplexVec Qhy; Qhy.reserve(NumTxAnts);
	ComplexMatrix2D invRh;
	ComplexVec HrxSymbs; HrxSymbs.reserve(NumTxAnts);
	double lambda = 0.000001;
	double lambdaML;
	Vecs pairlist = getNodeIdx(M, NumTxAnts);
	int vnum = pairlist.size();

	// cout<< "I'm before else!"<< endl;
	if (NumRxAnts >= NumTxAnts){  //lhq
		qrDecomposeEigen(H, Qh, R);

		// perform (Qy = Q'*rxSymbs); 
		for (auto i = 0; i < NumTxAnts; i++) {
			Qhy.push_back(inner_product(Qh[i].begin(), Qh[i].end(), rxSymbs.begin(), complex<double>(0.0)));
		} 
		lambdaML = SearchRadius * SearchRadius;
		}
	else {
		cholDecomposeEigen(H,invRh,Hh,R,lambda); // lhq  lambda = 0.000001
		// perform (Qy = inv(R') * H' * rxSymbs);
		for (auto i = 0; i < NumTxAnts; i++) {
			HrxSymbs.push_back(inner_product(Hh[i].begin(), Hh[i].end(), rxSymbs.begin(), complex<double>(0.0)));
		}
		for (auto i = 0; i < NumTxAnts; i++) {
			Qhy.push_back(inner_product(invRh[i].begin(), invRh[i].end(), HrxSymbs.begin(), complex<double>(0.0)));
		}
		lambdaML = SearchRadius * SearchRadius + lambda * NumTxAnts/M; // std::numeric_limits<double>::infinity();
		}


	// initalize node to root
	TreeNode node;
	// add the first symbol of the costellation for the last antenna
	node.level = NumTxAnts - 1;
	Constellation& cstltn = Cstltns[node.level];
	node.psv.push_front(cstltn.ConstlSymbs[0]);
	node.symIndices.push_front(0);
	node.bitLabels.insert(node.bitLabels.end(), cstltn.BitLabels[0].begin(), cstltn.BitLabels[0].end());
	for(int i = 0;i<NumTxAnts+1;i++){
		node.count1.push_front(0);
		node.count_1.push_front(0);
	} 
	node.count1[NumTxAnts - 1]=1;
	Vecd ped(NumTxAnts, 0.0); // partial euclidean distance
	double di; //distance increment

	// varibles related to tree traveral
	bool isTravelDone = false;
	NumNodesVisted = 0; // only count leaf nodes
	Matrix2Dd errorVector = {};
	int errorVectorNum = 0;

	while (!isTravelDone)
	{ 
		// if (TimeElapsed() >= cutoff_time)
        //     {
        //         return errorVector;
        //     }
		// compute distance increment
		complex<double> a = inner_product(node.psv.begin(), node.psv.end(), R[node.level].begin() + node.level, complex<double>(0.0));
		complex<double> b = Qhy[node.level] -a;
		di = norm(b);
		// cout<<"di="<<di<<endl;
		// compute partial Euclidean distance, none root node distance increamented di from its parent node.
		ped[node.level] = (node.level == NumTxAnts - 1) ? di : ped[node.level + 1] + di;
		Matrix2Di nodespair;
		Vecd errorDistance;
		bool notAllZero = false;
		
		if (node.level == 0) { // leaf node
			NumNodesVisted++;
			bool bML1 = ped.front() < lambdaML;
			bool bML2 = node.count_1[0]==node.count1[0];
			bool bML3 = node.count_1[0] <= M;
			if (bML1 && bML2 && bML3) { // smaller euclidean distance found
				
				bitLabelsML = Veci(node.bitLabels.begin(), node.bitLabels.end());
		
				for (int i:bitLabelsML){
					if (i != 0) {
						notAllZero =true;
						break;
					}
				}
					
				if (notAllZero){
					Matrix2Ds errorpair = getNodesPairIdx(bitLabelsML,M);
					double d = ped.front();
					
					for (auto i = 0; i < errorpair.size(); i++) {
						Vecd errorVector1 = {};
						errorVector1.push_back(d);
						// errorVector[errorVectorNum].push_back(d);
						for(string ii:errorpair[i]){
							int position = find(pairlist.begin(),pairlist.end(),ii)- pairlist.begin();
							if(position == vnum){
								cout<< "didn't find nodes!";
							}
							errorVector1.push_back(position+1);   // exist ???
						}
						errorVector.push_back(errorVector1);
					}
				}
			}
			isTravelDone = moveRight(node);
		}
		else 
		{	
			bool a1 = max(node.count1[node.level],node.count_1[node.level]) > M;
			bool a2 = abs(node.count1[node.level]-node.count_1[node.level]) > node.level;
			if (ped[node.level] >= lambdaML || a1 || a2)
				isTravelDone = moveRight(node); // trunc
			else
				moveDown(node);
		}
	
	}
	return errorVector;
}

void SphereDecoder::qrDecomposeEigen(const ComplexMatrix2D& H, ComplexMatrix2D& Qh, ComplexMatrix2D& R)
{
	//produces an economy-size decomposition, computes only the first N columns of Q and the first N rows of R.
	// Qh is hermitan tranpose of Q sized N-by-M;

	using namespace Eigen;

	Qh = ComplexMatrix2D(NumTxAnts); //hermitan tranpose of Q sized N-by-M
	R = ComplexMatrix2D(NumTxAnts);

	// Eigen domain perform QR decompostion
	typedef Matrix<complex<double>, Dynamic, Dynamic> cmat;
	cmat EigenH(NumRxAnts, NumTxAnts);
	for (auto i = 0; i < NumRxAnts; i++)
		for (auto j = 0; j < NumTxAnts; j++)
			EigenH(i, j) = H[i][j];

	HouseholderQR<cmat> qr(EigenH.rows(), EigenH.cols());
	qr.compute(EigenH);
	cmat EigenQ = qr.householderQ() * cmat::Identity(EigenH.rows(), EigenH.cols());
	cmat temp = qr.matrixQR().triangularView<Upper>();
	cmat EigenR = temp.topRows(EigenH.cols());

	// hermitan transpose
	cmat EigenQh = EigenQ.adjoint();

	//back to std namespace
	for (auto i = 0; i < NumTxAnts; i++) {
		Qh[i].reserve(NumRxAnts);
		for (auto j = 0; j < NumRxAnts; j++) {
			Qh[i].push_back(EigenQh(i, j));
		}
	}

	for (auto i = 0; i < NumTxAnts; i++) {
		R[i].reserve(NumTxAnts);
		for (auto j = 0; j < NumTxAnts; j++) {
			R[i].push_back(EigenR(i, j));
		}
	}
}

void SphereDecoder::cholDecomposeEigen(const ComplexMatrix2D& H, ComplexMatrix2D& invRh,ComplexMatrix2D& Hh, ComplexMatrix2D& R, double lambda = 0.000001)
{
	//produces an economy-size decomposition, computes only the first N columns of Q and the first N rows of R.
	// Qh is hermitan tranpose of Q sized N-by-M;

	using namespace Eigen;

	// Qh = ComplexMatrix2D(NumTxAnts); //hermitan tranpose of Q sized N-by-M
	invRh = ComplexMatrix2D(NumTxAnts); // inv(R') R is the upper matrix after cholesky factorization
	Hh = ComplexMatrix2D(NumTxAnts);
	R = ComplexMatrix2D(NumTxAnts);
	// Eigen domain perform QR decompostion
	typedef Matrix<complex<double>, Dynamic, Dynamic> cmat;
	cmat EigenH(NumRxAnts, NumTxAnts),EigenHh(NumTxAnts, NumRxAnts);
	for (int i = 0; i < NumRxAnts; i++){
		for (auto j = 0; j < NumTxAnts; j++){
			EigenH(i, j) = H[i][j];
		}
	}

	EigenHh = EigenH.adjoint();

	cmat ATA(NumTxAnts, NumTxAnts);
	ATA = EigenHh * EigenH + lambda * cmat::Identity(ATA.rows(), ATA.cols());

	LLT<cmat> lltOfA(ATA); // compute the Cholesky decomposition of A
	cmat EigenUpper = lltOfA.matrixU();  // R
	cmat EigenLower = EigenUpper.adjoint();//
	cmat invEigenLower = EigenLower.inverse();   // inv(R')

	//back to std namespace
	for (auto i = 0; i < NumTxAnts; i++) {
		Hh[i].reserve(NumTxAnts);
		for (auto j = 0; j < NumRxAnts; j++) {
			Hh[i].push_back(EigenHh(i, j));
		}
	}

	for (auto i = 0; i < NumTxAnts; i++) {
		R[i].reserve(NumTxAnts);
		for (auto j = 0; j < NumTxAnts; j++) {
			R[i].push_back(EigenUpper(i, j));
		}
	}

	
	
	for (auto i = 0; i < NumTxAnts; i++) {
		invRh[i].reserve(NumTxAnts);
		//cout << endl;
		for (auto j = 0; j < NumTxAnts; j++) {
			invRh[i].push_back(invEigenLower(i, j));
			// cout << " " << invEigenLower(i, j);
		}
	}
}

void SphereDecoder::moveDown(TreeNode& node)
{
	// node move done a level, always add the first symbol to the psv in the constellation of that level
	node.level--;
	Constellation& cstl = Cstltns[node.level];
	node.symIndices.push_front(0);
	node.psv.push_front(cstl.ConstlSymbs[0]);
	node.bitLabels.insert(node.bitLabels.begin(), cstl.BitLabels[0].begin(), cstl.BitLabels[0].end());
	node.count1[node.level] = node.count1[node.level+1] +1;
	node.count_1[node.level] = node.count_1[node.level+1];
};

bool SphereDecoder::moveRight(TreeNode& node)
{
	// try to move right, if already the rightest node of a subtree, then move up a level and go right
	while (true) {
		Constellation& cstl = Cstltns[node.level];
		bool a3 = ~(node.symIndices[0]==2 && node.count1[node.level+1]==0);
		if (node.symIndices[0] < cstl.M - 1 && a3) { // not yet to the far right,
			//level no change, the fist symbol of the psv will be changed to the next symbol in the constellation
			int symbIdx = node.symIndices[0] + 1;
			node.symIndices[0] = symbIdx;
			node.psv[0] = cstl.ConstlSymbs[symbIdx];
			node.bitLabels[0] = cstl.BitLabels[symbIdx][0];
			// for (auto i = 0; i < cstl.K; i++) node.bitLabels[i] = cstl.BitLabels[symbIdx][i];
			node.count_1[node.level]=node.count_1[node.level+1] + cstl.Count_1[node.symIndices[0]];
			node.count1[node.level]=node.count1[node.level+1] + cstl.Count1[node.symIndices[0]];
			return false;
		}
		else { // reach the rightest child of the parent node, UP
			if (node.level == NumTxAnts - 1) // root node
				return true; // traveral completed
			else { // go up a level
				node.psv.pop_front();
				node.symIndices.pop_front();
				node.bitLabels.pop_front();
				// for (auto i = 0; i < cstl.K; i++) node.bitLabels.pop_front();
				node.level++;
			}
		}
	}
};

Veci de2bi(int decNum, int NOut, string msb)
{
	// This function Convert decimal numbers to binary numbers.
	int tmp = decNum;
	Veci out(NOut, 0);
	if (msb == "right-msb") {
		int i = 0;
		while ((tmp > 0) & (i < NOut)) {
			out[i++] = tmp % 2;
			tmp /= 2;
		}
	}
	else {
		int i = NOut - 1;
		while ((tmp > 0) & (i >= 0)) {
			out[i--] = tmp % 2;
			tmp /= 2;
		}
	}

	return out;
}