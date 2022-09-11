#ifndef NR_SPHERE_DECODER
#define NR_SPHERE_DECODER

#include "nrModulation.h"
#include <string>
#include <complex>
#include <deque>
Vecs getNodeIdx(int M, int B);
Veci de2bi(int decNum, int NOut, std::string msb = "left-msb");
struct Constellation {
	std::string MoudType = "qpsk";
	int K = 2;
	int M = 4;
	Matrix2Di BitLabels{};
	ComplexVec ConstlSymbs{};
	Vecd Count1;
	Vecd Count_1;

	//constractor
	Constellation(std::string moduType) {
		std::transform(moduType.begin(), moduType.end(), moduType.begin(),
			[](unsigned char c) { return std::tolower(c); });
		MoudType = moduType;
		if (moduType == "bpsk") {
			M = 2; K = 1;
		}
		else if (moduType == "qpsk") {
			M = 4; K = 2;
		}
		else if (moduType == "3psk") {
			M = 3; 
		}
		else if (moduType == "16qam") {
			M = 16; K = 4;
		}
		else if (moduType == "64qam") {
			M = 64; K = 6;
		}
		else if (moduType == "256qam") {
			M = 256; K = 8;
		}
		else {
			M = 4; K = 2;
		}

		if (moduType == "3psk") {
			K = 3;
			BitLabels = {{1},{0},{-1}};
			ConstlSymbs = {{1},{0},{-1}};
			Count1 = {{1},{0},{0}};
			Count_1 = {{0},{0},{1}};
		}
		else{
			BitLabels = Matrix2Di(M);
			ConstlSymbs = ComplexVec(M);
			for (auto i = 0; i < M; i++) {
				BitLabels[i] = de2bi(i, K);
				ConstlSymbs[i] = nrModuMapper(BitLabels[i], moduType, M)[0];
			}
		}
	}
};
struct TreeNode {
	int level = 0; // 0 ~ Nt-1, Nt is number of Tx antennas, root level is Nt-1, leaf nodes with level 0
	std::deque<std::complex<double> > psv{}; // each node has a corresponding partial symbol vector
	std::deque<int> symIndices{}; // the symbIndices corrsponding to each symbol in the psv
	std::deque<int> bitLabels{}; // the bitlabels corrsponding to each symbol in the psv
	std::deque<int> count1{};
	std::deque<int> count_1{};
};

class SphereDecoder {
private:
	int NumRxAnts = 0; // number of rx antennas
	int NumTxAnts = 0;	// number of transmit antennas
	// constellations for each tx antenna, as when 2 codewords are supported ,
	// different antennas can have different modulations
	std::vector<Constellation> Cstltns{};
	Veci Ks; // records number of bits for each decoded tx symbols on each tx antenna
	int NumTotBits = 0;; // total number of bits.
	Veci Offsets{}; // records the starting position for bit labels for each decoded tx symbols
	int NumNodesVisted = 0; // number of nodes visted in sphere decoding
	double SearchRadius = 0.2; // Search radius to control tree pruning level
public:
	SphereDecoder(int nRxAnts, int nTxAnts, std::vector<std::string>& moduTypes, double searchRadius = 0.2);
	~SphereDecoder();
	Vecd operator()(ComplexMatrix2D& H, ComplexVec& rxSymbs); // soft decoding
	Matrix2Dd hardSphereDecode(ComplexMatrix2D& H, ComplexVec& rxSymbs, int M, double SearchRadius=0.9);
	static Matrix2Ds getNodesPairIdx(Veci E, int M);

public: //geters
	int getNumOfNodeVisited() { return NumNodesVisted; };
	Veci getKs() { return Ks; };
	Veci getOffsets() { return Offsets; };
	int getTotNumBits() { return NumTotBits; }

private:
	void nodeInit(TreeNode& node);
	void qrDecomposeEigen(const ComplexMatrix2D& H, ComplexMatrix2D& Qh, ComplexMatrix2D& R);
	void cholDecomposeEigen(const ComplexMatrix2D& H, ComplexMatrix2D& invRh,ComplexMatrix2D& Hh,ComplexMatrix2D& R,double lambda);
	void moveDown(TreeNode& node);
	bool moveRight(TreeNode& node);
};

#endif // !NR_SPHERE_DECODER
