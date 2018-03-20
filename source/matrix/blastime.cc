/* Ergo, version 3.6, a program for linear scaling electronic structure
 * calculations.
 * Copyright (C) 2017 Elias Rudberg, Emanuel H. Rubensson, Pawel Salek,
 * and Anastasia Kruchinina.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Primary academic reference:
 * Kohn−Sham Density Functional Theory Electronic Structure Calculations 
 * with Linearly Scaling Computational Time and Memory Usage,
 * Elias Rudberg, Emanuel H. Rubensson, and Pawel Salek,
 * J. Chem. Theory Comput. 7, 340 (2011),
 * <http://dx.doi.org/10.1021/ct100611z>
 * 
 * For further information about Ergo, see <http://www.ergoscf.org>.
 */

/** @file blastime.cc

    \brief Benchmark code for matrix-matrix multiplication (gemm)
    operation, measuring timings for different matrix sizes.
*/

#include <iostream>
#include <fstream>
#include <iomanip> /* For setprecision in fstream */
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <sys/time.h>
#include <vector>

#include "mat_gblas.h"

static const int MIN_TIME_PER_STEP = 5;
static const int SIZE_INCREMENT = 2;

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

template<class T>
static void tomatlabfile(char* name,T* values,int s,std::ofstream& output);

template<typename real>
int mainFun(int maxDim, double* timev, double* gflops, bool writeTomFile) {
  try {
    //    bool CPUtime = false;
    const real ONE=1.0;
    const real ZERO=0.0;
    double start, end;
    int i;
    int steps = 500000;//50000000;

    /* Find reasonable number of steps */
    double secondsTaken = 0;
    int testSize = SIZE_INCREMENT;
    while (secondsTaken < MIN_TIME_PER_STEP) {
      steps = steps*3;
      real* A=new real [testSize*testSize];
      real* B=new real [testSize*testSize];
      real* C=new real [testSize*testSize];
      for(i = testSize*testSize-1; i>=0; i--) A[i] = 1.0;
      for(i = testSize*testSize-1; i>=0; i--) B[i] = 1.0;
      for(i = testSize*testSize-1; i>=0; i--) C[i] = 0.0;
      int m = testSize;
      int n = testSize;
      int k = testSize;
      start = get_wall_seconds();
      for(int j=0; j<steps; j++) {
	       mat::gemm("N","N",&m,&n,&k,&ONE,A,&m,B,&k,&ZERO,C,&m);
      }
      end = get_wall_seconds();
      secondsTaken = end-start;
    }
    printf("%d tests took %6.2f seconds.\n", steps, secondsTaken);
    
    /* Run actual benchmark */
    int maxStep = maxDim/SIZE_INCREMENT;
    for (int step = 1; step <= maxStep; step++) {
      int size = step*SIZE_INCREMENT;
      real* A=new real [size*size];
      real* B=new real [size*size];
      real* C=new real [size*size];
      for(i = size*size-1; i>=0; i--) A[i] = 1.0;
      for(i = size*size-1; i>=0; i--) B[i] = 1.0;
      for(i = size*size-1; i>=0; i--) C[i] = 0.0;
      int m = size;
      int n = size;
      int k = size;
      
      start = get_wall_seconds();
      for(int j=0; j<steps; j++) {
	       mat::gemm("N","N",&m,&n,&k,&ONE,A,&m,B,&k,&ZERO,C,&m);
      }
      end = get_wall_seconds();
      secondsTaken = end-start;
      timev[step-1] = secondsTaken;
      gflops[step-1]=(2*pow(double(size),3)+4*pow(double(size),2))/(timev[step-1]*1e9)*steps;

      // gflops[step-1]=((2*pow(size,3)+4*pow(size,2))/timev[step-1])*(steps/1000000000);
      delete[] A;
      delete[] B;
      delete[] C;
      std::cout<<"size="<<std::setw(4)<<size
	       <<"    steps="<<std::setw(8)<<steps
	       <<"    time="<<std::setw(6)<<timev[step-1]
	       <<"    Gflops="<<std::setw(10)
	       <<gflops[step-1]<<std::endl;
      if (timev[step-1]>MIN_TIME_PER_STEP*2) {
	/* This prediction does not really work for large relative
	   matrix size increments. But it does its job sufficiently well. */
	int newSteps = 
	  ceil(double(steps)*MIN_TIME_PER_STEP*1.5/double(timev[step-1]));
	std::cout << "Recomputing new steps "<< newSteps <<" = "
		  << steps << " * " << MIN_TIME_PER_STEP << " / "
		  << timev[step-1] << std::endl;
	steps = newSteps;
      }
    }
  }
  catch (std::exception e) {
    std::cout << "Exception caught: "<<e.what() << std::endl;
    std::exit(1);
  }
  return 0;
}


template<typename T>
static void tomatlabfile(const char* name, T* values, int s,
			 std::ofstream& output) {
  output<<name<<"=[";
  for (int i=0;i<s;i++)
    {
      output<<std::setprecision(10)<<values[i]<<'\n';
    }
  output<<"];"<<std::endl;
}




int main(int argc,char* argv[]) {
  int maxDim;
  char path[200];
  bool writeTomFile = true;
  switch (argc) {
  case 2:
    maxDim = atoi(argv[1]); /* Max matrix dimension. */
    writeTomFile = false;
    break;
  case 3:
    maxDim = atoi(argv[1]); /* Max matrix dimension. */
    strcpy(path, argv[2]);  /* Matlab filename.      */
    break;
  default:
    std::cerr<<"Wrong number of input arguments. Give at least one argument (max matrix dimension) and optionally as a second argument the name of the m-file to create." << std::endl;
    std::exit(1);
  }
  int maxSize = maxDim/SIZE_INCREMENT;
  std::vector<double> timevDouble(maxSize);
  std::vector<double> gflopsDouble(maxSize);
  std::vector<double> timevSingle(maxSize);
  std::vector<double> gflopsSingle(maxSize);
  maxDim = maxSize*SIZE_INCREMENT;

  time_t startTime;
  time_t endTime;
  std::cout<<"Starting gemm benchmark, double precision"<<std::endl;
  time(&startTime);
  if (!mainFun<double>(maxDim, &timevDouble[0], &gflopsDouble[0], writeTomFile)) {
    time(&endTime);
    std::cout<<"Ended gemm benchmark, double precision, wall time: "
	     <<endTime-startTime
	     <<std::endl;    
  }
  std::cout<<"Starting gemm benchmark, single precision"<<std::endl;
  time(&startTime);
  if (!mainFun<float>(maxDim, &timevSingle[0], &gflopsSingle[0], writeTomFile)) {
    time(&endTime);
    std::cout<<"Ended gemm benchmark, single precision, wall time: "
	     <<endTime-startTime
	     <<std::endl;
  }
  
  if (writeTomFile) {
    std::ofstream output(path);
    if (!output) {
      std::cout<<"Cannot open outputfile"<<std::endl;
      std::exit(1);
    }
    output<<"nv=1:" << SIZE_INCREMENT <<":"<<maxDim<<";"<<std::endl;
    tomatlabfile<double>("timeDouble"  , &timevDouble [0], maxSize, output);
    tomatlabfile<double>("GflopsDouble", &gflopsDouble[0], maxSize, output);
    tomatlabfile<double>("timeSingle"  , &timevSingle [0], maxSize, output);
    tomatlabfile<double>("GflopsSingle", &gflopsSingle[0], maxSize, output);
    output<<"minX = 0;\n"
	  <<"maxX = max(nv);\n"
	  <<"minY = 0;\n"
	  <<"maxY = max([GflopsDouble GflopsSingle]);\n"
	  <<"lwidth = 1;\n"
	  <<"fsize  = 12;\n"
	  <<"figure;\n"
	  <<std::endl
	  <<"subplot(221)\n"
	  <<"plot(nv,GflopsDouble,'LineWidth',lwidth);nflops=100000;\n"
	  <<"axis([minX maxX minY maxY(1)])\n"
	  <<"set(gca,'FontSize',fsize)\n"
	  <<"title('Double Precision')\n"
	  <<"xlabel('Matrix Size')\n"
	  <<"ylabel('Gflops')\n"
	  <<"grid on\n"
	  <<std::endl
	  <<"subplot(222)\n"
	  <<"plot(nv,GflopsSingle,'LineWidth',lwidth);nflops=100000;\n"
	  <<"axis([minX maxX minY maxY(2)])\n"
	  <<"set(gca,'FontSize',fsize)\n"
	  <<"title('Single Precision')\n"
	  <<"xlabel('Matrix Size')\n"
	  <<"ylabel('Gflops')\n"
	  <<"grid on\n";
  }
  
  
  std::exit(0);
}
