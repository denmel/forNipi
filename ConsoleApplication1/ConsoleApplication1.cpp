#include <iostream>
#include <fstream>
#include <math.h>
#include "Matrix.h"
#include "Calc.h"
#include <vector>

using namespace std;
const string inputFile = "C:\\Users\\denis\\source\\repos\\ConsoleApplication1\\x64\\Debug\\input.txt";
const string outputFile = "C:\\Users\\denis\\source\\repos\\ConsoleApplication1\\x64\\Debug\\output.txt";
int main(int argc, char* argv[])
{	
	Calc calc = Calc(0.1, 0.000001);
	calc.readData(inputFile);
	calc.save(outputFile);
	while (calc.t < 20) {	
		//calc.stepEuler();
		calc.stepRK();
		//calc.stepSymplecticVerle();
		//calc.stepImplicitEuler();	
		//calc.stepSymplecticEuler();
		
		calc.save(outputFile);
		calc.print();
	}
}