#include "utils.h"

void Utils::printMatrix(int **pMat, int rows, int cols)
{
	int i, j;
	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			cout << pMat[i][j];
		}
		cout << endl;
	}
	cout << endl;
}

int Utils::gf2rank(int **A, int m, int n)
{
	int csum = 0;
	int i, j, k, l, tmp;
	int r = 0;
	int **tmpA = new int*[m];
	for (i = 0; i < m; i++)
	{
		tmpA[i] = new int[n];
		memcpy(tmpA[i], A[i], sizeof(int)*n);
	}
	for (i = 0; i < m - 1; i++)
	{
		if (tmpA[i][i] == 0)
		{
			for (j = i + 1; j < m; j++)
			{
				if (tmpA[j][i] != 0)
				{
					for (k = 0; k < n; k++)
					{
						tmpA[i][k] ^= tmpA[j][k];
					}
					break;
				}
			}
		}
		if (tmpA[i][i] == 0) continue;
		for (j = i + 1; j < m; j++)
		{
			if (tmpA[j][i] == 1)
			{
				for (k = i; k < n; k++)
				{
					tmpA[j][k] = tmpA[j][k] ^ tmpA[i][k];
				}
			}
		}
	}
	for (k = n - 1; k >= 1; k--)
	{
		for (i = k; i >= 1; i--)
		{
			if (tmpA[i][k] == 1) break;
		}
		if (i == 0) continue;
		for (j = i - 1; j >= 0; j--)
		{
			if (tmpA[j][k] == 1)
			{
				for (l = i; l <= k; l++)
				{
					tmpA[j][l] = tmpA[j][l] ^ tmpA[i][l];
				}
			}
		}
	}
	for (i = 0; i < m; i++)
	{
		csum = 0;
		for (j = 0; j < n; j++)
		{
			csum += tmpA[i][j];
		}
		if (csum != 0) r++;
	}
	for (i = 0; i < m; i++)
	{
		delete[] tmpA[i];
	}
	delete[] tmpA;
	return r;
}

void Utils::blockXor(unsigned long long *dest, unsigned long long *src, int longwords)
{
	/*char *end = src + bytes;
	unsigned long long *pStart = (unsigned long long*)src;
    unsigned long long *pEnd = (unsigned long long*)end;
	unsigned long long *pDest = (unsigned long long*)dest;
	while (pStart < pEnd) {
		*pDest = ((*pDest) ^ (*pStart));
		pStart++;
		pDest++;
	}*/

	for (int i = 0; i < longwords; i++)
	{
		dest[i] = (dest[i] ^ src[i]);
	}
}

int Utils::getFileLength(const char *pFileName)
{
	ifstream fin(pFileName, ios::binary | ios::in);
	fin.seekg(0, ios::beg);
	streamoff start = fin.tellg();
	fin.seekg(0, ios::end);
	streamoff end = fin.tellg();
	int filesize = end - start;
	fin.close();
	return filesize;
}

void Utils::IntArrayXor(int *arr1, int *arr2, int num)
{
	int *end = arr2 + num;
	unsigned __int64 *pStart = (unsigned __int64*)arr2;
	unsigned __int64 *pEnd = (unsigned __int64*)end;
	unsigned __int64 *pDest = (unsigned __int64*)arr1;
	while (pStart < pEnd) {
		*pDest = ((*pDest) ^ (*pStart));
		pStart++;
		pDest++;
	}
}

void Utils::inverse(int **A, int **invA, int dim)
{
	int i, j, k, tmp;
	int **B = new int*[dim];
	for (i = 0; i < dim; i++){
		B[i] = new int[2 * dim];
		memset(B[i], 0, sizeof(int)* 2 * dim);
		memcpy(B[i], A[i], sizeof(int)*dim);
		B[i][i + dim] = 1;
	}
	for (i = 0; i < dim - 1; i++)
	{
		if (B[i][i] == 0)
		{
			for (j = i + 1; j < dim; j++)
			{
				if (B[j][i] != 0)
				{
					for (k = 0; k < 2 * dim; k++)
					{
						tmp = B[i][k];
						B[i][k] = B[j][k];
						B[j][k] = tmp;
					}
					break;
				}
			}
		}
		for (j = i + 1; j < dim; j++)
		{
			if (B[j][i] == 1)
			{
				for (k = i; k < 2 * dim; k++)
				{
					B[j][k] = B[j][k] ^ B[i][k];
				}
			}
		}
	}
	for (i = dim - 1; i > 0; i--){
		for (j = i - 1; j >= 0; j--){
			if (B[j][i] == 1){
				for (k = 0; k < 2 * dim; k++){
					B[j][k] = B[j][k] ^ B[i][k];
				}
			}
		}
	}
	for (i = 0; i < dim; i++){
		memcpy(invA[i], *(B + i) + dim, sizeof(int)*dim);
	}
	for (i = 0; i < dim; i++){
		delete[] B[i];
		B[i] = NULL;
	}
	delete[]B;
	B = NULL;
}

void Utils::fileCompare(const char *file1, int erasedInd)
{
	string sFile = string(file1);
	int dotPos = sFile.find('.');
	string sFilePrefix = string(file1, file1 + dotPos);
	string sFileSuffix = string(file1 + dotPos);
	ostringstream oss;
	oss << erasedInd;
	ifstream fin;
	string sfile = "E:\\RegeneratingCodeExperiment\\RegeneratingCodeExperiment\\Coding\\" + sFilePrefix
		            + "_b" + oss.str() + "recovered" + sFileSuffix;
	int fileLen = getFileLength(sfile.c_str());
	fin.open(sfile.c_str(), ios::in | ios::binary);
	char *buf1 = new char[fileLen + 1];
	fin.read(buf1, fileLen);
	buf1[fileLen] = '\0';
	fin.close();
	sfile = "E:\\RegeneratingCodeExperiment\\RegeneratingCodeExperiment\\backup\\" + sFilePrefix
		+ "_b" + oss.str() + sFileSuffix;
	fin.open(sfile.c_str(), ios::in | ios::binary);
	char *buf2 = new char[fileLen + 1];
	fin.read(buf2, fileLen);
	buf2[fileLen] = '\0';
	fin.close();
	if (!memcmp(buf1, buf2, fileLen))
	{
		cout << string(sFilePrefix + "_b" + oss.str() + sFileSuffix) << "和"
			 << string(sFilePrefix + "_b" + oss.str() + "recovered" + sFileSuffix) << "相同." << endl;
	}
	else
	{
		cout << string(sFilePrefix + "_b" + oss.str() + sFileSuffix) << "和"
			 << string(sFilePrefix + "_b" + oss.str() + "recovered" + sFileSuffix) << "不相同." << endl;
	}
	delete[] buf1;
	delete[] buf2;
}