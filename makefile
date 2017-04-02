Rand_Mat_Reg_Code:Rand_Mat_Reg_Code.o utils.o
	g++ -o Rand_Mat_Reg_Code Rand_Mat_Reg_Code.o utils.o

Rand_Mat_Reg_Code.o:Rand_Mat_Reg_Code.cpp Rand_Mat_Reg_Code.h
	g++ -c Rand_Mat_Reg_Code.cpp

utils.o:utils.cpp utils.h
	g++ -c utils.cpp

clean:
	rm utils.o Rand_Mat_Reg_Code.o Rand_Mat_Reg_Code
