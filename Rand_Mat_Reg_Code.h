#ifndef _Rand_Mat_Reg_Code_H_
#define _Rand_Mat_Reg_Code_H_
#include <iostream>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <ctime>
#include <string>
#include <sstream>
#include <set>
#include <unistd.h>
using namespace std;

class RegeneratingCode
{
public:
	RegeneratingCode(int n, int k, int d, int l, string category);

	void GenerateCodingMatrix(const char *pFile);

	void EncodeStripFile(unsigned long long **coding, int offset, int rowSize);

	void ReadSourceFile(const char* pFileName, char *buffer, int filesize);

	void MbrStripCopyFile(int *offset);

	void MsrStripCopyFile(int *offset);

	void Encode(const char *pFileName);

	void WriteChunksToFile(const char *pFileName, unsigned long long **coding, int chunkSize);

	void RandomEraseBlock();

	void GenerateHelpNodeIndex(const char *pCodingMatFile);

	void RegenerateChunkFile(const char *pFileName, char **helpData, int readIns);

	void RegenerateChunkFile(const char *pFileName, int **invMatrix, unsigned long long **helpData, int readIns, int cols);

	void Regenerate(const char *pFileName, const char *pCodingMatFile);

	~RegeneratingCode();

	int erasedInd;
//private:
	int                n, k, d, l;       //再生码编码参数
	int                **codingMatrix;
	int                **recoverMatrix;
	int                **invRecoverMatrix;
	unsigned long long **data;
	char               *fileBuffer;
	string             category;
	set<int>           helpNodeIndexs;
};
#endif
