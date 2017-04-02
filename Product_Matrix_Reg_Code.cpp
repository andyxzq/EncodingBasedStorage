/********************************************
Date:2016/11/30
Author:xuzhiqiang
Decription:Experiment about regenerating code using sparse random matrix vandermonde matrix
********************************************/
#include "Product_Matrix_Reg_Code.h"

RegeneratingCode::RegeneratingCode(int n, int k, int d, int delta, string category)
{
	this->n = n;
	this->k = k;
	this->d = d;
	this->l = d + delta;
	this->category = category;
	cout << "********************************" << endl;
	cout << "Regenerating Code's Parameters (n, k, d, l) as follows:" << endl;
	cout << "n = " << n << endl;
	cout << "k = " << k << endl;
	cout << "d = " << d << endl;
	cout << "l = " << (d+delta) << endl;
	cout << "code category = " << category << endl;
	cout << "********************************" << endl;
	int i;
	recoverMatrix = new int*[d];
	for(i = 0; i < d; i++)
	{
		recoverMatrix[i] = new int[d];
	}
	if(category == "SR-MBR" || category == "SR-MSR")
	{
		codingMatrix = new int*[n];
		
		invRecoverMatrix = new int*[d];
		for (i = 0; i < n; i++)
		{
			if (i < d) 
			{
				invRecoverMatrix[i] = new int[d];
				memset(invRecoverMatrix[i], 0, d);
			}
			codingMatrix[i] = new int[d];
		}
	}
	else
	{
		vandermondeMatrix = new int[n * d];
		if(category == "Van-MSR")
		{
			subMatrix = new int[n * d / 2];
		    diagonalMatrix = new int[n * n];
	    }
	    recoverArray = new int[d * d];
	}
	if (category == "SR-MBR")
	{
		data1 = new unsigned long long*[d];
		for (i = 0; i < d; i++)
		{
			data1[i] = new unsigned long long[d];
			memset(data1[i], 0, d*8);
		}
	}
	else if (category == "SR-MSR")
	{
		data1 = new unsigned long long*[2*(d - k + 1)];
		for (i = 0; i < 2*(d-k+1); i++)
		{
			data1[i] = new unsigned long long[d-k+1];
			memset(data1[i], 0, (d - k + 1)*8);
		}
	}
	else if(category == "Van-MBR")
	{
		data2 = new char*[d];
		for (i = 0; i < d; i++)
		{
			data2[i] = new char[d];
			memset(data2[i], 0, d);
		}
		parity = new char*[n - d];
		for (i = 0; i < n - d; i++)
		{
			parity[i] = new char[d];
			memset(parity[i], 0, d);
		}
	}
	else
	{
		data2 = new char*[2*(d - k + 1)];
		for (i = 0; i < 2 * (d - k + 1); i++)
		{
			data2[i] = new char[d - k + 1];
			memset(data2[i], 0, d - k + 1);
		}
		parity = new char*[n - d];
		for (i = 0; i < n - d; i++)
		{
			parity[i] = new char[d-k+1];
			memset(parity[i], 0, d-k+1);
		}
	}
	helpNodeIndexs.clear();
}

//assign the ith row's weight as lambda
void RegeneratingCode::AssignRow(int **A, int i, int dim, int lambda, int *positions)
{
	int pos, flag, j, l;
	memset(positions, 0, sizeof(int)*lambda);
	for (j = 0; j < lambda; j++)
	{
		pos = rand() % dim;
		positions[j] = pos;
		if (j == 0) A[i][pos] = 1;
		else
		{
			flag = 0;
			while (!flag)
			{
				for (l = 0; l < j; l++)
				{
					if (positions[l] == positions[j])
					{
						break;
					}
				}
				if (l == j)
				{
					flag = 1;
					A[i][pos] = 1;
				}
				else
				{
					pos = rand() % dim;
					positions[j] = pos;
				}
			}
		}
	}
}

int Equals(int* A, int* B, int len)
{
	int i;
	for (i = 0; i < len; i++)
	{
		if (A[i] != B[i]) return 0;
	}
	return 1;
}

void RegeneratingCode::GenerateCodingMatrix(const char *pFile)
{
	if (n <= d)
	{
		cerr << "n always greater than d!" << endl;
		exit(1);
	}
	int i, j, flag;
	char *curDir = new char[200];
	getcwd(curDir, 200);
	mkdir("Coding", S_IRWXU|S_IRWXG|S_IROTH);
	string destFile = string(curDir) + string("/Coding/") + string(pFile);
	ofstream fout(destFile.c_str());

    int rowweight = 2 * (int)(log10(d) / log10(2)) + 1;    //row weight of coding matrix
    int *positions = new int[rowweight];
	if(category == "SR-MBR" || category == "SR-MSR")
	{
		flag = 0;
		srand((unsigned int)time(NULL));
		while (!flag)
		{
			/*for (i = 0; i < n; i++)
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
			}*/
			for (i = 0; i < n; i++)
			{
				memset(codingMatrix[i], 0, sizeof(int)*d);
				if (i < d) codingMatrix[i][i] = 1;
			}
			for (i = d; i < n; i++)
			{
				AssignRow(codingMatrix, i, d, rowweight, positions);
				if (i > d)
				{
					flag = 0;
					while (!flag)
					{
						for (j = 0; j < i; j++)
						{
							if (Equals(codingMatrix[i], codingMatrix[j], d) == 1)
							{
								flag = 1;
								break;
							}
						}
						if (flag == 1)
						{
							memset(codingMatrix[i], 0, sizeof(int)*d);
							AssignRow(codingMatrix, i, d, rowweight, positions);
							flag = 0;
						}
						else if (flag == 0)
						{
							flag = 1;
						}
					}
				}
			}
			if (Utils::gf2rank(codingMatrix, n, d) == d) flag = 1;
		}
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < d; j++)
			{
				fout << codingMatrix[i][j] << " ";
			}
			fout << "\n";
		}
	}
	else
	{
        vandermondeMatrix = reed_sol_big_vandermonde_distribution_matrix(n, d, 8);
        if(category == "Van-MSR")
        {
	        subMatrix = reed_sol_big_vandermonde_distribution_matrix(n, d / 2, 8);
	        memset(diagonalMatrix, 0, sizeof(int)*n*n);
	        i = 0;
	        srand((unsigned int)time(NULL));
	        flag = 0;
	        while(i < n)
	        {
	        	diagonalMatrix[i * (n + 1)] = n - i;//rand() % 255 + 1;
	        	/*for(j = 0; j < i; j++)
	        	{
	        		if(diagonalMatrix[j * (n + 1)] == diagonalMatrix[i * (n + 1)])
	        		{
	        			flag = 1;
	        			break;
	        		}
	        	}
	        	if(flag == 0)
	        	{
	        		i++;
	        	}
	        	else
	        	{
	        		flag = 0;
	        	}*/
	        	i++;
	        }
	        int *product = jerasure_matrix_multiply(diagonalMatrix, subMatrix, n, n, n, d / 2, 8);
	        for(i = 0; i < n; i++)
	        {
	        	memcpy(vandermondeMatrix + i * d, subMatrix + i * d / 2, sizeof(int) * d / 2);
	        	memcpy(vandermondeMatrix + i * d + d / 2, product + i * d / 2, sizeof(int) * d / 2);
	        }
	    }
        for(i = 0; i < n * d; i++)
        {
        	if((i != 0) && (i % d == 0)) fout<<"\n";
        	fout<<vandermondeMatrix[i]<<" ";
        }
	}
	fout.close();
}

//encode the file strip
void RegeneratingCode::EncodeStripFile(unsigned long long **coding, int offset, int rowSize)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		if (i < d)
		{
			memcpy(coding[i]+offset, data1[i], rowSize*8);
		}
		else
		{
			memset(coding[i]+offset, 0, rowSize*8);
			for (j = 0; j < d; j++)
			{
				if (codingMatrix[i][j] == 1)
				{
					Utils::blockXor(coding[i]+offset, data1[j], rowSize);
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

void RegeneratingCode::MbrStripCopyFile(int *offset)
{
	int i, j;
	for (i = 0; i < k; i++)
	{
		if(category == "SR-MBR")
		{
		    memcpy(data1[i], fileBuffer + *offset, (i + 1)*8);
	    }
	    else if(category == "Van-MBR")
	    {
	    	memcpy(data2[i], fileBuffer + *offset, i + 1);
	    }
		*offset += i + 1;
	}
	if (d > k)
	{
		for (i = k; i < d; i++)
		{
			if(category == "SR-MBR")
			{
			    memcpy(data1[i], fileBuffer + *offset, k*8);
			}
			else if(category == "Van-MBR")
			{
				memcpy(data2[i], fileBuffer + *offset, k);
			}
			*offset += k;
		}
	}
	for (i = 0; i < d; i++)
	{
		for (j = i + 1; j < d; j++)
		{
			if(category == "SR-MBR")
			{
			    data1[i][j] = data1[j][i];
			}
			else if(category == "Van-MBR")
			{
				data2[i][j] = data2[j][i];
			}
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
			if(category == "SR-MSR")
			{			
				memcpy(data1[i], fileBuffer + *offset, (i + 1)*8);
			}
			else
			{
                memcpy(data2[i], fileBuffer + *offset, i + 1);
			}
			*offset += i + 1;
		}
		else
		{
			if(category == "SR-MSR")
			{
			    memcpy(data1[i], fileBuffer + *offset, (i - d/2 + 1)*8);
		    }
		    else
		    {
		    	memcpy(data2[i], fileBuffer + *offset, i - d/2 + 1);
		    }
			*offset += i - d / 2 + 1;
		}
	}
	for (i = 0; i < d; i++)
	{
		if (i < d / 2)
		{
			for (j = i + 1; j < d / 2; j++)
			{
				if(category == "SR-MSR")
				{
				    data1[i][j] = data1[j][i];
				}
				else
				{
					data2[i][j] = data2[j][i];
				}
			}
		}
		else
		{
			for (j = i - d / 2 + 1; j < d / 2; j++)
			{
				if(category == "SR-MSR")
				{
				    data1[i][j] = data1[j + d / 2][i - d/2];
			    }
			    else
			    {
                    data2[i][j] = data2[j + d / 2][i - d/2];
			    }
			}
		}
	}
}

//wtite SR-MBR or SR-MSR coded block into disk files
void RegeneratingCode::WriteChunksToFile1(const char *pFileName, unsigned long long **coding, int chunkSize)
{
	string sFile = string(pFileName);
	int dotPos = sFile.find('.');
	string sFilePrefix = string(pFileName, pFileName + dotPos);
	string sFileSuffix = string(pFileName + dotPos);
	char *curDir = new char[200];
	getcwd(curDir, 200);
	string destDir = string(curDir) + string("/Coding/");
	for (int i = 0; i < n; i++)
	{
		ostringstream oss;
		oss << i+1;
		string sFileName = destDir + sFilePrefix + "_b" + oss.str() + sFileSuffix;
		ofstream fout(sFileName.c_str(), ios::binary|ios::out|ios::app);
		fout.write((char*)coding[i], chunkSize*8);
		fout.close();
	}
	delete[] curDir;
	curDir = NULL;
}

//wtite Van-MBR or Van-MSR coded block into disk files
void RegeneratingCode::WriteChunksToFile2(const char *pFileName, char **coding, int chunkSize)
{
	string sFile = string(pFileName);
	int dotPos = sFile.find('.');
	string sFilePrefix = string(pFileName, pFileName + dotPos);
	string sFileSuffix = string(pFileName + dotPos);
	char *curDir = new char[200];
	getcwd(curDir, 200);
	string destDir = string(curDir) + string("/Coding/");
	for (int i = 0; i < n; i++)
	{
		ostringstream oss;
		oss << i+1;
		string sFileName = destDir + sFilePrefix + "_b" + oss.str() + sFileSuffix;
		ofstream fout(sFileName.c_str(), ios::binary|ios::out|ios::app);
		fout.write(coding[i], chunkSize);
		fout.close();
	}
	delete[] curDir;
	curDir = NULL;
}

void RegeneratingCode::Encode(const char *pFileName)
{
	int i, j, readIns;
	int filesize = Utils::getFileLength(pFileName);
	int stripSize;                                                               
	if (category == "SR-MBR" || category == "Van-MBR")
	{
		stripSize = (2 * k * d - k * k + k) / 2;                              
	}
	else if (category == "SR-MSR" || category == "Van-MSR")
	{
		stripSize = (d - k + 1)*(d - k + 2);
	}
	//cout << "Original File Size: " << filesize <<"B"<< endl;
	//cout << "File Strip Size: " << stripSize * 8<<"B"<< endl;
	if(category == "SR-MBR" || category == "SR-MSR")
	{
	    readIns = filesize / (stripSize * 8) + ((filesize % (stripSize * 8) == 0) ? 0 : 1);
    }
    else
    {
    	readIns = filesize / stripSize + ((filesize % stripSize == 0) ? 0 : 1);
    }
	//cout << "Operation Repeated Times: " << readIns << endl;
    
    if(category == "SR-MBR" || category == "SR-MSR")
    {
	    fileBuffer = new char[readIns * stripSize * 8];                   
	    memset(fileBuffer, 0, readIns * stripSize * 8);
	}
	else
	{
		fileBuffer = new char[readIns * stripSize];                   
	    memset(fileBuffer, 0, readIns * stripSize);
	}
	ReadSourceFile(pFileName, fileBuffer, filesize);                     

	int alpha;                                                                  
	unsigned long long **coding;                        //store coded block of SR-MBR or SR-MSR
	char **bufferedCoding;                              //store coded block of Van-MBR or Van-MSR                 
	if (category == "SR-MBR" || category == "Van-MBR") alpha = d;
	else if(category == "SR-MSR"|| category == "Van-MSR") alpha = d - k + 1;
	if (category == "SR-MBR" || category == "SR-MSR")
	{
		coding = new unsigned long long*[n];                                  
		for (i = 0; i < n; i++)
		{
			coding[i] = new unsigned long long[alpha*readIns + 1];
			memset(coding[i], 0, (alpha*readIns + 1)*8);
		}
	}
	else
	{
		bufferedCoding = new char*[n];
        for(i = 0; i < n; i++)
        {
	        bufferedCoding[i] = new char[alpha*readIns + 1];
	        memset(bufferedCoding[i], 0, alpha*readIns+1);
        }
	}

	int cnt = 1;
	int offset = 0;                                                             
	cout << "Encoding Process Start..." << endl;
	clock_t startTime = clock();
	while (cnt <= readIns)
	{
		if (category == "SR-MBR")
		{
			for (i = 0; i < d; i++)
			{
				memset(data1[i], 0, d*8);
			}
			MbrStripCopyFile(&offset);
		}
		else if(category == "Van-MBR")
		{
			for (i = 0; i < d; i++)
			{
				memset(data2[i], 0, d);
			}
			MbrStripCopyFile(&offset);
		}
		else if(category == "SR-MSR")
		{
			for (i = 0; i < 2 * (d - k + 1); i++)
			{
				memset(data1[i], 0, (d - k + 1)*8);
			}
			MsrStripCopyFile(&offset);
		}
		else
		{
			for (i = 0; i < 2 * (d - k + 1); i++)
			{
				memset(data2[i], 0, d - k + 1);
			}
			MsrStripCopyFile(&offset);
		}

		if(category == "SR-MBR" || category == "SR-MSR")
		{
		    EncodeStripFile(coding, (cnt - 1)*alpha, alpha);                     
		}
		else
		{
            jerasure_matrix_encode(d, n - d, 8, vandermondeMatrix, data2, parity, alpha);
            for(i = 0; i < n; i++)
            {
            	if(i < d)
            	{
            		memcpy(bufferedCoding[i] + cnt * alpha, data2[i], alpha);
            	}
            	else
            	{
            		memcpy(bufferedCoding[i] + cnt * alpha, parity[i - d], alpha);
            	}
            }

		}
		++cnt;
	}

	if(category == "SR-MBR" || category == "SR-MSR")
	{
	    WriteChunksToFile1(pFileName, coding, alpha*readIns);                      
	}
	else
	{
		WriteChunksToFile2(pFileName, bufferedCoding, alpha*readIns);
	}
	clock_t endTime = clock();
	cout << "Encoding Process End!" << endl;
	cout << "Time: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
}

void RegeneratingCode::GenerateHelpNodeIndex(const char *pCodingMatFile)
{
	int i, j;
	ifstream fin(pCodingMatFile);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < d; j++)
		{
			if(category == "SR-MBR" || category == "SR-MSR")
			{
			    fin >> codingMatrix[i][j];
			}
			else
			{
				fin >> vandermondeMatrix[i * d + j];
			}
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
			if(category == "SR-MBR" || category == "SR-MSR")
			{
			    memcpy(recoverMatrix[i++], codingMatrix[*iter - 1], sizeof(int)*d);
		    }
		    else
		    {
		    	memcpy(recoverMatrix[i], vandermondeMatrix+(*iter-1)*d, sizeof(int)*d);
		    	memcpy(recoverArray + i*d, vandermondeMatrix+(*iter-1)*d, sizeof(int)*d);
		    	i++;
		    }
		}
		if(category == "SR-MSR" || category == "SR-MBR")
		{
			if (Utils::gf2rank(recoverMatrix, d, d) == d) 
			{
				flag = 1;
			}
			else helpNodeIndexs.clear();
			cout << ".";
	    }
	    else if(category == "Van-MSR")
	    {
	    	if(jerasure_invertible_matrix(recoverArray, d, 8) == 1)
	    	{
	    		flag = 1;
	    	}
	    	else helpNodeIndexs.clear();
	    	cout<<".";
	    }
	    else break;
	}
	cout << endl;
}

void RegeneratingCode::RegenerateChunkFile(const char *pFileName, int **invMatrix, unsigned long long **helpData, int readIns, int cols)
{
	int i, j, t;
	unsigned long long **regeneratedData = new unsigned long long*[readIns];
	for (i = 0; i < readIns; i++)
	{
		regeneratedData[i] = new unsigned long long[cols];
		memset(regeneratedData[i], 0, cols*8);
	}
	for (i = 0; i < readIns; i++)
	{
		for (j = 0; j < d; j++)
		{
			for (t = 0; t < d; t++)
			{
				if (invMatrix[j][t] == 1)
				{
					regeneratedData[i][j] = regeneratedData[i][j] ^ helpData[i][t];
				}
			}
		}
		for (j = d; j < cols; j++)
		{
			for (t = d; t < cols; t++)
			{
				if (invMatrix[j - d][t - d] == 1)
				{
					regeneratedData[i][j] = regeneratedData[i][j] ^ helpData[i][t];
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
	getcwd(curDir, 200);
	string targetFile = string(curDir) + string("/Coding/") + sFilePrefix + "_b" + oss.str() +
		"recovered" + sFileSuffix;
	ofstream fout(targetFile.c_str(), ios::app | ios::binary | ios::out);
	for (i = 0; i < readIns; i++)
	{
		if (category == "SR-MBR")
		{
			fout.write((char*)regeneratedData[i], d * 8);
		}
		else if (category == "SR-MSR")
		{
			Utils::blockXor(regeneratedData[i], regeneratedData[i] + 3 * d / 2, d / 2);
			fout.write((char*)regeneratedData[i], d * 4);
		}
	}
	fout.close();
}

void RegeneratingCode::Regenerate(const char *pFileName, const char *pCodingMatFile)
{
	int i, j, t, offset, m;
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
	char *curDir = new char[200];
	getcwd(curDir, 200);
	string targetFile = string(curDir) + string("/Coding/") + sFilePrefix + "_b" + oss.str() + sFileSuffix;
	int fileChunkSize = Utils::getFileLength(targetFile.c_str());
	fin.open(targetFile.c_str(), ios::in | ios::binary);
	char *tmpBuffer = new char[fileChunkSize + 1];
	fin.read(tmpBuffer, fileChunkSize);
	fin.close();
	string destFile = string(curDir) + string("/backup/") + sFilePrefix + "_b" + oss.str() + sFileSuffix;
	fout.open(destFile.c_str(), ios::out | ios::binary);
	fout.write(tmpBuffer, fileChunkSize);
	fout.close();
	remove(targetFile.c_str());                                              
	cout << "Deleted Node Id: " << "b" + oss.str() << endl;

	string sCodingMatFile = string(curDir) + string("/Coding/") + string(pCodingMatFile);
	GenerateHelpNodeIndex(sCodingMatFile.c_str());
	int alpha, readIns;
	if (category == "SR-MBR" || category == "Van-MBR") alpha = d;
	else if(category == "SR-MSR" || category == "Van-MSR") alpha = d - k + 1;
	if(category == "SR-MBR" || category == "SR-MSR")
	{
	    readIns = (fileChunkSize % (alpha*8) == 0 ? fileChunkSize / (alpha*8) : fileChunkSize / (alpha*8) + 1);
    }
    else
    {
    	readIns = (fileChunkSize % alpha == 0 ? fileChunkSize / alpha : fileChunkSize / alpha + 1);
    }
	cout << "Operations Repeated Times: " << readIns << endl;
	if (category == "SR-MBR")
	{
		cout << "Repair Bandwidth: " << l*readIns*8 << "B" << endl;
	}
	else if (category == "SR-MSR")
	{
		cout << "Repair Bandwidth: " << 2*l*readIns*8 << "B" <<endl;
	}
	else
	{
        cout << "Repair Bandwidth: " << d*readIns << "B" << endl;
	}

	unsigned long long **coding = new unsigned long long*[d];                              
	char **bufferedCoding = new char*[d];
	for (i = 0; i < d; i++)
	{
		if(category == "SR-MBR" || category == "SR-MSR")
		{
		    coding[i] = new unsigned long long[alpha*readIns];
		    memset(coding[i], 0, alpha*readIns*8);
	    }
	    else
	    {
	    	bufferedCoding[i] = new char[alpha*readIns];
	    	memset(bufferedCoding[i], 0, alpha*readIns);
	    }
	}
	cout<<"1.++++++++++"<<endl;

	i = 0;
	for (iter = helpNodeIndexs.begin(); iter != helpNodeIndexs.end(); ++iter)
	{
		oss.str("");
		oss << *iter;
		targetFile = string(curDir) + string("/Coding/") + sFilePrefix + "_b" + oss.str() + sFileSuffix;
		fin.open(targetFile.c_str(), ios::in | ios::binary);
		if(category == "SR-MBR" || category == "SR-MSR")
		{
		    fin.read((char*)coding[i++], alpha*readIns*8);
	    }
	    else
	    {
	    	fin.read(bufferedCoding[i++], alpha*readIns);
	    }
		fin.close();
	}

	unsigned long long **helpData1 = new unsigned long long*[readIns];
	char **helpData2 = new char*[d];
	int *invRecoverArray, **decodingMatrix;                      
	int cols;
	if (category == "SR-MBR" || category == "Van-MBR" || category == "Van-MSR")
	{
		cols = d;
	}
	else if (category == "SR-MSR")
	{
		cols = 2 * d;
	}
	if(category == "SR-MBR" || category == "SR-MSR")
	{
		for (i = 0; i < readIns; i++)
		{
		    helpData1[i] = new unsigned long long[cols];
	        memset(helpData1[i], 0, cols*8);
		}
        Utils::inverse(recoverMatrix, invRecoverMatrix, d);
        targetFile = string(curDir) + string("/Coding/") + "DecodingMatrix.txt";
        fout.open(targetFile.c_str(), ios::out);
        for(i = 0; i < d; i++)
        {
        	for(j = 0; j < d; j++)
        	{
        		fout<<invRecoverMatrix[i][j];
        	}
        	fout<<"\n";
        }
        fout.close();
    }
    else
    {
    	invRecoverArray = new int[d*d];
    	decodingMatrix = new int*[d];
        for(i = 0; i < d; i++)
        {
        	helpData2[i] = new char[alpha];
        	memcpy(recoverArray+i*d, recoverMatrix[i], d*sizeof(int));
        }
        jerasure_invert_matrix(recoverArray, invRecoverArray, d, 8);
        targetFile = string(curDir) + string("/Coding/") + "DecodingMatrix.txt";
        fout.open(targetFile.c_str(), ios::out);
        for(i = 0; i < d; i++)
        {
        	decodingMatrix[i] = new int[d];
        	memcpy(decodingMatrix[i], invRecoverArray + i*d, d*sizeof(int));
        	for(j = 0; j < d; j++)
        	{
        		fout<<invRecoverArray[i + j * d];
        	}
        	fout<<"\n";
        }
        fout.close();
    }
    
    unsigned long long **regeneratedData = new unsigned long long*[readIns];
	for (i = 0; i < readIns; i++)
	{
		regeneratedData[i] = new unsigned long long[cols];
		memset(regeneratedData[i], 0, cols*8);
	}

	cout << "Regenerating Start..." << endl;
	clock_t startTime = clock();
	if(category == "SR-MBR" || category == "SR-MSR")
	{
		offset = 0;
		//help nodes encode chunks locally
		for (i = 0; i < readIns; i++)
		{
			for (j = 0; j < d; j++)
			{
				for (t = 0; t < alpha; t++)                                 
				{
					if (codingMatrix[erasedInd - 1][t] == 1)
					{
						helpData1[i][j] = helpData1[i][j] ^ coding[j][offset + t];
					}
				}
				
			}
			for (j = d; j < cols; j++)
			{
				for (t = alpha; t < d; t++)
				{
					if (codingMatrix[erasedInd - 1][t] == 1)
					{
						helpData1[i][j] = helpData1[i][j] ^ coding[j - d][offset + t - alpha];
					}
				}
			}
			offset += alpha;
		}
		//RegenerateChunkFile(pFileName, invRecoverMatrix, helpData1, readIns, cols);
		//Replace node regenerate lost chunk using received data
		for (i = 0; i < readIns; i++)
		{
			for (j = 0; j < d; j++)
			{
				for (t = 0; t < d; t++)
				{
					if (invRecoverMatrix[j][t] == 1)
					{
						regeneratedData[i][j] = regeneratedData[i][j] ^ helpData1[i][t];
					}
				}
			}
			for (j = d; j < cols; j++)
			{
				for (t = d; t < cols; t++)
				{
					if (invRecoverMatrix[j - d][t - d] == 1)
					{
						regeneratedData[i][j] = regeneratedData[i][j] ^ helpData1[i][t];
					}
				}
			}
		}
		
		targetFile = string(curDir) + string("/Coding/") + sFilePrefix + "_b" + oss.str() +
			"recovered" + sFileSuffix;
		fout.open(targetFile.c_str(), ios::app | ios::binary | ios::out);
		for (i = 0; i < readIns; i++)
		{
			if (category == "SR-MBR")
			{
				fout.write((char*)regeneratedData[i], d * 8);
			}
			else if (category == "SR-MSR")
			{
				Utils::blockXor(regeneratedData[i], regeneratedData[i] + 3 * d / 2, d / 2);
				fout.write((char*)regeneratedData[i], d * 4);
			}
		}
		fout.close();
	}
	else
	{
		targetFile = string(curDir) + string("/Coding/") + sFilePrefix + "_b" + oss.str() +
			"recovered" + sFileSuffix;
		fout.open(targetFile.c_str(), ios::app | ios::binary | ios::out);
		offset = 0;
		int *matrix1 = new int[(d+1)*d];
		int *matrix2 = new int[2*d*d];
		memset(matrix1, 0, (d+1)*d*sizeof(int));
		memset(matrix2, 0, 2*d*d*sizeof(int));
		char** strip1 = new char*[1];
		strip1[0] = new char[d];
		char** strip2 = new char*[d];
		char** strip3 = new char*[d];
		char *strip4 = new char[alpha + 1];   //regenerate strip data of Van-MSR
		for(i = 0; i < d; i++)
		{
			strip2[i] = new char[1];
			strip3[i] = new char[1];
			matrix1[i*d + i] = 1;	
			matrix2[i*d + i] = 1;
		}
		for(i = d; i < 2*d; i++)
		{
			memcpy(matrix2+i*d, decodingMatrix[i-d], d*sizeof(int));
		}
        memcpy(matrix1 + d*d, vandermondeMatrix + (erasedInd - 1)*d, d*sizeof(int));
		for(i = 0; i < readIns; i++)
		{
			for(j = 0; j < d; j++)
			{
				memset(helpData2[j], 0, alpha);
				memcpy(helpData2[j], bufferedCoding[j]+offset, alpha);
			}
            jerasure_matrix_encode(d, 1, 8, matrix1, helpData2, strip1, alpha);
            for(m = 0; m < d; m++)
            {
            	strip2[m][0] = strip1[0][m];
            }
            jerasure_matrix_encode(d, d, 8, matrix2, strip2, strip3, 1);
            for(m = 0; m < d; m++)
            {
            	strip1[0][m] = strip3[m][0];
            }
            if(category == "Van-MSR")
            {
            	char *product = new char[alpha];
            	galois_w08_region_multiply(strip1[0] + alpha, diagonalMatrix[(n+1)*(erasedInd-1)], alpha, product, 0);
            	galois_region_xor(product, strip1[0], alpha);
            }
            if(category == "Van-MBR")
            {
                fout.write(strip1[0], d);
            }
            else
            {
            	fout.write(strip1[0], alpha);
            }
			offset += alpha;
		}
		fout.close();
	}
	clock_t endTime = clock();
	cout << "Lost Chunk Regenerated Success!" << endl;
	cout << "Time: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	for(i = 0; i < readIns; i++) delete[] regeneratedData[i];
	delete[] regeneratedData;
}

RegeneratingCode::~RegeneratingCode()
{
	int i;
	for (i = 0; i < d; i++)
	{
		delete[] recoverMatrix[i];
	}
	delete[] recoverMatrix;
	delete[] fileBuffer;
}

int main()
{
	const char* srcFile = "140MB.mp4";
	const char* codingMatFile = "CodingMatrix.txt";
	RegeneratingCode rc(100, 40, 78, 20, "Van-MSR");
	rc.GenerateCodingMatrix(codingMatFile);
	rc.Encode(srcFile);
	rc.Regenerate(srcFile, codingMatFile);
	Utils::fileCompare(srcFile, rc.erasedInd);
	return 0;
}