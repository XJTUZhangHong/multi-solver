#include "output.h"

void output1d(Fluid1d* fluids, Block1d block)
{
	ofstream ResultFile;
	int N = block.nodex, ghost = block.ghost;
	const char* filePath_U = "../data/U.txt";
	const char* filePath_R = "../data/R.txt";

	ResultFile.open(filePath_U);
	for (int i = ghost; i < N + ghost; i++)
	{
		ResultFile << (fluids[i].comp1.convar[1] + fluids[i].comp2.convar[1]) / (fluids[i].comp1.convar[0] + fluids[i].comp2.convar[0]) << endl;
	}
	ResultFile.close();

	ResultFile.open(filePath_R);
	for (int i = ghost; i < N + ghost; i++)
	{
		ResultFile << (fluids[i].comp1.convar[0] + fluids[i].comp2.convar[0]) << endl;
	}
	ResultFile.close();
}