#ifndef _UTILS_H_
#define _UTILS_H_
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

class Utils
{
public:
	static void printMatrix(int **pMat, int rows, int cols);
	static int gf2rank(int **A, int m, int n);
	static void blockXor(unsigned long long *dest, unsigned long long *src, int longwords);
	static int getFileLength(const char *pFileName);
	static void fileCompare(const char *file1, int erasedInd);
	static void IntArrayXor(int *arr1, int *arr2, int num);
	static void inverse(int **A, int **invA, int dim);
};
#endif