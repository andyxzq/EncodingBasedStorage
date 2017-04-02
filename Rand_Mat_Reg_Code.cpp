/********************************************
Date:2016/11/30
Author:xuzhiqiang
Decription:Experiment about sparse random matrix regenerating code(SR-MBR and SR-MSR)
********************************************/
#include "Rand_Mat_Reg_Code.h"
#include "utils.h"

RegeneratingCode::RegeneratingCode(int n, int k, int d, int delta, string category)
{
	this->n = n;
	this->k = k;
	this->d = d;
	this->l = d + delta;
	this->category = category;
	cout << "********************************" << endl;
	cout << "再生码参数(n, k, d, l)设定如下:" << endl;
	cout << "n = " << n << endl;
	cout << "k = " << k << endl;
	cout << "d = " << d << endl;
	cout << "l = " << (d+delta) << endl;
	cout << "code category = " << category << endl;
	cout << "********************************" << endl;
	int i;
	codingMatrix = new int*[n];
	recoverMatrix = new int*[d];
	invRecoverMatrix = new int*[d];
	for (i = 0; i < n; i++)
	{
		if (i < d) 
		{
			recoverMatrix[i] = new int[d];
			invRecoverMatrix[i] = new int[d];
			memset(invRecoverMatrix[i], 0, d);
		}
		codingMatrix[i] = new int[d];
	}
	if (category == "SR-MBR")
	{
		data = new unsigned long long*[d];
		for (i = 0; i < d; i++)
		{
			data[i] = new unsigned long long[d];
			memset(data[i], 0, d*8);
		}
	}
	else if (category == "SR-MSR")
	{
		data = new unsigned long long*[2*(d - k + 1)];
		for (i = 0; i < 2*(d-k+1); i++)
		{
			data[i] = new unsigned long long[d-k+1];
			memset(data[i], 0, (d - k + 1)*8);
		}
	}
	helpNodeIndexs.clear();
}

/*生成编码矩阵*/
void RegeneratingCode::GenerateCodingMatrix(const char *pFile)
{
	if (n <= d)
	{
		cerr << "参数n应该大于参数d" << endl;
		exit(1);
	}
	int i, j;
	char *curDir = new char[100];
	_getcwd(curDir, 100);
	cout << curDir << endl;
	string destFile = string(curDir) + string("\\Coding\\") + string(pFile);
	ofstream fout(destFile);
	int flag = 0;
	srand(time(NULL));
	while (!flag)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < d; j++)
			{
				if (i < d)
				{
					if (j == i) codingMatrix[i][j] = 1;
					else codingMatrix[i][j] = 0;
				}
				else
				{
					codingMatrix[i][j] = rand() % 2;
				}
			}
		}
		if (Utils::gf2rank(codingMatrix, n, d) == d) flag = 1;
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < d; j++)
		{
			if (j != 0) fout << " ";
			fout << codingMatrix[i][j];
		}
		fout << "\n";
	}
	fout.close();
}

/*对数据块条带进行编码*/
void RegeneratingCode::EncodeStripFile(unsigned long long **coding, int offset, int rowSize)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		if (i < d)
		{
			//memcpy(coding[i]+offset, data[i], rowSize*8);
			memcpy(coding[i], data[i], rowSize * 8);
		}
		else
		{
			//memset(coding[i]+offset, 0, rowSize*8);
			memset(coding[i], 0, rowSize * 8);
			for (j = 0; j < d; j++)
			{
				if (codingMatrix[i][j] == 1)
				{
					//Utils::blockXor(coding[i]+offset, data[j], rowSize);
					Utils::blockXor(coding[i], data[j], rowSize);
				}
			}
		}
	}
}

void RegeneratingCode::ReadSourceFile(const char* pFileName, char *buffer, int filesize)
{
	memset(buffer, 0, filesize);
	ifstream fin(pFileName, ios::binary | ios::in);
	fin.read(buffer, filesize);
	fin.close();
}

/*拷贝数据至数据块矩阵*/
void RegeneratingCode::MbrStripCopyFile(int *offset)
{
	int i, j;
	for (i = 0; i < k; i++)
	{
		memcpy(data[i], fileBuffer + *offset, (i + 1)*8);
		*offset += i + 1;
	}
	if (d > k)
	{
		for (i = k; i < d; i++)
		{
			memcpy(data[i], fileBuffer + *offset, k*8);
			*offset += k;
		}
	}
	for (i = 0; i < d; i++)
	{
		for (j = i + 1; j < d; j++)
		{
			data[i][j] = data[j][i];
		}
	}
}

void RegeneratingCode::MsrStripCopyFile(int *offset)
{
	int i, j;
	for (i = 0; i < d; i++)
	{
		if (i < d / 2)
		{
			memcpy(data[i], fileBuffer + *offset, (i + 1)*8);
			*offset += i + 1;
		}
		else
		{
			memcpy(data[i], fileBuffer + *offset, (i - d/2 + 1)*8);
			*offset += i - d / 2 + 1;
		}
	}
	for (i = 0; i < d; i++)
	{
		if (i < d / 2)
		{
			for (j = i + 1; j < d / 2; j++)
			{
				data[i][j] = data[j][i];
			}
		}
		else
		{
			for (j = i - d / 2 + 1; j < d / 2; j++)
			{
				data[i][j] = data[j + d / 2][i - d/2];
			}
		}
	}
}

/*写编码块到文件*/
void RegeneratingCode::WriteChunksToFile(const char *pFileName, unsigned long long **coding, int chunkSize)
{
	string sFile = string(pFileName);
	int dotPos = sFile.find('.');
	string sFilePrefix = string(pFileName, pFileName + dotPos);
	string sFileSuffix = string(pFileName + dotPos);
	char *curDir = new char[100];
	_getcwd(curDir, 100);
	string destDir = string(curDir) + string("\\Coding\\");
	for (int i = 0; i < n; i++)
	{
		ostringstream oss;
		oss << i+1;
		string sFileName = destDir + sFilePrefix + "_b" + oss.str() + sFileSuffix;
		ofstream fout(sFileName, ios::binary|ios::out|ios::app);
		fout.write((char*)coding[i], chunkSize*8);
		fout.close();
	}
	delete[] curDir;
	curDir = NULL;
}

void RegeneratingCode::Encode(const char *pFileName)
{
	int i, j;
	int filesize = Utils::getFileLength(pFileName);
	int stripSize;                                                               //单次读取的数据条块大小
	if (category == "SR-MBR")
	{
		stripSize = (2 * k * d - k * k + k) / 2;                              
	}
	else if (category == "SR-MSR")
	{
		stripSize = (d - k + 1)*(d - k + 2);
	}
	cout << "原文件大小: " << filesize <<"B"<< endl;
	cout << "数据条块大小: " << stripSize * 8<<"B"<< endl;
	int readIns = filesize / (stripSize * 8) + ((filesize % (stripSize * 8) == 0) ? 0 : 1); //记录总的读取次数
	cout << "操作重复次数: " << readIns << endl;

	fileBuffer = new char[readIns * stripSize * 8];                             //文件缓存
	memset(fileBuffer, 0, readIns * stripSize * 8);
	ReadSourceFile(pFileName, fileBuffer, filesize);                            //读取文件至内存缓存区

	unsigned long long **coding = new unsigned long long*[n];                   //存储文件编码块
	int alpha;                                                                  //单次生成的编码条块大小
	if (category == "SR-MBR") alpha = d;
	else if(category == "SR-MSR") alpha = d - k + 1;
	for (i = 0; i < n; i++)
	{
		//coding[i] = new unsigned long long[alpha*readIns + 1];
		//memset(coding[i], 0, (alpha*readIns + 1)*8);
		coding[i] = new unsigned long long[alpha + 1];
		memset(coding[i], 0, (alpha + 1) * 8);
	}

	int cnt = 1;
	int offset = 0;                                                             //记录每次从buffer开始读取的位置
	cout << "编码过程开始..." << endl;
	clock_t startTime = clock();
	while (cnt <= readIns)
	{
		if (category == "MBR")
		{
			for (i = 0; i < d; i++)
			{
				memset(data[i], 0, d*8);
			}
			MbrStripCopyFile(&offset);
		}
		else
		{
			for (i = 0; i < 2 * (d - k + 1); i++)
			{
				memset(data[i], 0, (d - k + 1)*8);
			}
			MsrStripCopyFile(&offset);
		}
		EncodeStripFile(coding, (cnt - 1)*alpha, alpha);                         //对每个strip文件进行编码
		WriteChunksToFile(pFileName, coding, alpha);
		++cnt;
	}
	//WriteChunksToFile(pFileName, coding, alpha*readIns);                         //将编码后的块写到磁盘文件
	clock_t endTime = clock();
	cout << "编码过程结束!" << endl;
	cout << "Time: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
//	for (i = 0; i < n; i++) delete[] coding[i];
//	delete[] coding;
}

void RegeneratingCode::RegenerateChunkFile(const char *pFileName, char **helpData, int readIns)
{
	int i, j, q, t, u;
	int **tmpMat = new int*[l];
	for (i = 0; i < l; i++)
	{
		tmpMat[i] = new int[d];
		memcpy(tmpMat[i], recoverMatrix[i], sizeof(int)*d);
	}
	for (u = 0; u < readIns; u++)
	{
		for (i = 0; i < l - 1; i++)
		{
			if (recoverMatrix[i][i] == 0)
			{
				for (j = i + 1; j < l; j++)
				{
					if (recoverMatrix[j][i] != 0)
					{
						Utils::IntArrayXor(recoverMatrix[i], recoverMatrix[j], d);
						helpData[u][i] = helpData[u][i] ^ helpData[u][j];
						break;
					}
				}
			}
			if (recoverMatrix[i][i] == 0) continue;
			for (j = i + 1; j < l; j++)
			{
				if (recoverMatrix[j][i] == 1)
				{
					Utils::IntArrayXor(recoverMatrix[j], recoverMatrix[i], d);
					helpData[u][j] = helpData[u][j] ^ helpData[u][i];
				}
			}
		}
		for (q = d - 1; q >= 1; q--)
		{
			for (i = q; i >= 1; i--)
			{
				if (recoverMatrix[i][q] == 1) break;
			}
			if (i == 0) continue;
			for (j = i - 1; j >= 0; j--)
			{
				if (recoverMatrix[j][q] == 1)
				{
					for (t = i; t <= q; t++)
					{
						recoverMatrix[j][t] ^= recoverMatrix[i][t];
					}
					helpData[u][j] = helpData[u][j] ^ helpData[u][i];
				}
			}
		}
		for (i = 0; i < l; i++)
		{
			memcpy(recoverMatrix[i], tmpMat[i], sizeof(int)*d);
		}
	}
	string sFile = string(pFileName);
	int dotPos = sFile.find('.');
	string sFilePrefix = string(pFileName, pFileName + dotPos);
	string sFileSuffix = string(pFileName + dotPos);
	ostringstream oss;
	oss << erasedInd;
	char *curDir = new char[200];
	_getcwd(curDir, 200);
	string targetFile = string(curDir) + string("\\Coding\\") + sFilePrefix + "_b" + oss.str() +
		                "recovered" + sFileSuffix;
	ofstream fout(targetFile, ios::app | ios::binary | ios::out);
	for (i = 0; i < readIns; i++)
	{
		fout.write(helpData[i], d);
	}
	fout.close();
	for (i = 0; i < l; i++)
	{
		delete[] tmpMat[i];
	}
	delete[] tmpMat;
	delete[] curDir;
}

/*生成帮助恢复节点的索引*/
void RegeneratingCode::GenerateHelpNodeIndex(const char *pCodingMatFile)
{
	int i, j;
	ifstream fin(pCodingMatFile);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < d; j++)
		{
			fin >> codingMatrix[i][j];
		}
	}
	fin.close();

	set<int>::iterator iter;
	int flag = 0;
	while (!flag)
	{
		i = 0;
		srand(time(NULL));
		while (i < d)
		{
			int index = rand() % n + 1;
			if (helpNodeIndexs.find(index) == helpNodeIndexs.end() && index != erasedInd)
			{
				helpNodeIndexs.insert(index);
				++i;
			}
		}
		i = 0;
		for (iter = helpNodeIndexs.begin(); iter != helpNodeIndexs.end(); iter++)
		{
			memcpy(recoverMatrix[i++], codingMatrix[*iter - 1], sizeof(int)*d);
		}
		if (Utils::gf2rank(recoverMatrix, d, d) == d) 
		{
			flag = 1;
		}
		else helpNodeIndexs.clear();
		cout << ".";
	}
	cout << endl;
}

/*恢复丢失的数据块*/
void RegeneratingCode::RegenerateChunkFile(const char *pFileName, int **invMatrix, unsigned long long **helpData, int readIns, int cols)
{
	int i, j, t;
	unsigned long long **regenerateData = new unsigned long long*[readIns];
	for (i = 0; i < readIns; i++)
	{
		regenerateData[i] = new unsigned long long[cols];
		memset(regenerateData[i], 0, cols*8);
	}
	for (i = 0; i < readIns; i++)
	{
		for (j = 0; j < d; j++)
		{
			for (t = 0; t < d; t++)
			{
				if (invMatrix[j][t] == 1)
				{
					regenerateData[i][j] = regenerateData[i][j] ^ helpData[i][t];
				}
			}
		}
		for (j = d; j < cols; j++)
		{
			for (t = d; t < cols; t++)
			{
				if (invMatrix[j - d][t - d] == 1)
				{
					regenerateData[i][j] = regenerateData[i][j] ^ helpData[i][t];
				}
			}
		}
	}
	string sFile = string(pFileName);
	int dotPos = sFile.find('.');
	string sFilePrefix = string(pFileName, pFileName + dotPos);
	string sFileSuffix = string(pFileName + dotPos);
	ostringstream oss;
	oss << erasedInd;
	char *curDir = new char[200];
	_getcwd(curDir, 200);
	string targetFile = string(curDir) + string("\\Coding\\") + sFilePrefix + "_b" + oss.str() +
		"recovered" + sFileSuffix;
	ofstream fout(targetFile, ios::app | ios::binary | ios::out);
	for (i = 0; i < readIns; i++)
	{
		if (category == "SR-MBR")
		{
			fout.write((char*)regenerateData[i], d*8);
		}
		else if (category == "SR-MSR")
		{
			Utils::blockXor(regenerateData[i], regenerateData[i] + 3 * d / 2, d / 2);
			fout.write((char*)regenerateData[i], d * 8 / 2);
		}
	}
	fout.close();
	//int fileChunkSize = Utils::getFileLength(targetFile.c_str());
	//cout << "file chunk size: " << fileChunkSize << endl;
}

void RegeneratingCode::Regenerate(const char *pFileName, const char *pCodingMatFile)
{
	int i, j, t, offset;
	ifstream fin;
	ofstream fout;
	set<int>::iterator iter;
	srand(time(NULL));
	erasedInd = rand()%n + 1;

	string sFile = string(pFileName);
	int dotPos = sFile.find('.');
	string sFilePrefix = string(pFileName, pFileName + dotPos);
	string sFileSuffix = string(pFileName + dotPos);
	ostringstream oss;
	oss << erasedInd;
	char *curDir = new char[100];
	_getcwd(curDir, 100);
	string targetFile = string(curDir) + string("\\Coding\\") + sFilePrefix + "_b" + oss.str() + sFileSuffix;
	int fileChunkSize = Utils::getFileLength(targetFile.c_str());
	//cout << "file chunk size: " << fileChunkSize << endl;
	fin.open(targetFile, ios::in | ios::binary);
	char *tmpBuffer = new char[fileChunkSize + 1];
	fin.read(tmpBuffer, fileChunkSize);
	fin.close();
	string destFile = string(curDir) + string("\\backup\\") + sFilePrefix + "_b" + oss.str() + sFileSuffix;
	fout.open(destFile, ios::out | ios::binary);
	fout.write(tmpBuffer, fileChunkSize);
	fout.close();
	remove(targetFile.c_str());                                              //删除随机选定的文件块
	cout << "删除节点编号为: " << "b" + oss.str() << endl;

	string sCodingMatFile = string(curDir) + string("\\Coding\\") + string(pCodingMatFile);
	GenerateHelpNodeIndex(sCodingMatFile.c_str());
	int alpha;
	if (category == "SR-MBR") alpha = d;
	else if(category == "SR-MSR") alpha = d - k + 1;
	int readIns = (fileChunkSize % (alpha*8) == 0 ? fileChunkSize / (alpha*8) : fileChunkSize / (alpha*8) + 1);
	cout << "操作重复次数为: " << readIns << endl;
	if (category == "SR-MBR")
	{
		cout << "修复带宽: " << l*readIns*8 << "B" << endl;
	}
	else if (category == "SR-MSR")
	{
		cout << "修复带宽: " << 2*l*readIns*8 << "B" <<endl;
	}

	unsigned long long **coding = new unsigned long long*[d];                              //存储帮助恢复的文件块
	for (i = 0; i < d; i++)
	{
		coding[i] = new unsigned long long[alpha*readIns];
		memset(coding[i], 0, alpha*readIns*8);
	}
	i = 0;
	for (iter = helpNodeIndexs.begin(); iter != helpNodeIndexs.end(); ++iter)
	{
		oss.str("");
		oss << *iter;
		targetFile = string(curDir) + string("\\Coding\\") + sFilePrefix + "_b" + oss.str() + sFileSuffix;
		fin.open(targetFile, ios::in | ios::binary);
		fin.read((char*)coding[i++], alpha*readIns*8);
		fin.close();
	}
	unsigned long long **helpData = new unsigned long long*[readIns];                       //存储帮助节点需要发送的数据
	int cols;
	if (category == "SR-MBR")
	{
		cols = d;
	}
	else if (category == "SR-MSR")
	{
		cols = 2 * d;
	}
	for (i = 0; i < readIns; i++)
	{
		helpData[i] = new unsigned long long[cols];
		memset(helpData[i], 0, cols*8);
	}

	Utils::inverse(recoverMatrix, invRecoverMatrix, d);
	cout << "开始恢复丢失块..." << endl;
	clock_t startTime = clock();
	offset = 0;
	for (i = 0; i < readIns; i++)
	{
		for (j = 0; j < d; j++)
		{
			for (t = 0; t < alpha; t++)                                 
			{
				if (codingMatrix[erasedInd - 1][t] == 1)
				{
					helpData[i][j] = helpData[i][j] ^ coding[j][offset + t];
				}
			}
			
		}
		for (j = d; j < cols; j++)
		{
			for (t = alpha; t < d; t++)
			{
				if (codingMatrix[erasedInd - 1][t] == 1)
				{
					helpData[i][j] = helpData[i][j] ^ coding[j - d][offset + t - alpha];
				}
			}
		}
		offset += alpha;
	}
	RegenerateChunkFile(pFileName, invRecoverMatrix, helpData, readIns, cols);
	clock_t endTime = clock();
	cout << "丢失的数据块回复成功!" << endl;
	cout << "Time: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	for (i = 0; i < d; i++)
	{
		delete[] coding[i];
	}
	for (i = 0; i < readIns; i++)
	{
		delete[] helpData[i];
	}
	delete[] helpData;
	delete[] coding;
	delete[] tmpBuffer;
}

RegeneratingCode::~RegeneratingCode()
{
	int i;
	for (i = 0; i < n; i++)
	{
		if (i < d) delete[] data[i];
		if (i < d) 
		{
			delete[] recoverMatrix[i];
			delete[] invRecoverMatrix[i];
		}
		delete[] codingMatrix[i];
	}
	delete[] codingMatrix;
	delete[] recoverMatrix;
	delete[] invRecoverMatrix;
	delete[] data;
}

int main()
{
	const char* srcFile = "140MB.mp4";
	const char* codingMatFile = "CodingMatrix.txt";
	RegeneratingCode rc(100, 50, 60, 75, "SR-MBR");
	rc.GenerateCodingMatrix(codingMatFile);
	rc.Encode(srcFile);
	rc.Regenerate(srcFile, codingMatFile);
	Utils::fileCompare(srcFile, rc.erasedInd);
	return 0;
}
