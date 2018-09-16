#include <iostream>
#include <random>
#include <omp.h>

#define MATSIZEX 128
#define MATSIZEY 128
#define RANGE_HIGH 31
#define RANGE_LOW 0

void populate_matrix(unsigned char matrix[MATSIZEX][MATSIZEY])
{
    std::random_device rd;
    std::mt19937 gen(rd());//mesenne twister engine
    std::uniform_int_distribution<unsigned char> dis(RANGE_LOW, RANGE_HIGH);
	
	#pragma omp parallel
	#pragma omp for
    for (int ii=0; ii<MATSIZEX; ++ii)
	{
		for(int jj=0; jj< MATSIZEY; ++jj)
		{
			matrix[ii][jj] = dis(gen);
		}
	}
}



void inline count( unsigned char matrix[MATSIZEX][MATSIZEY], size_t count_vals[] )
{

	#pragma omp parallel
	{
	#pragma omp for
	for(int kk=RANGE_LOW; kk<RANGE_HIGH+1; kk++)
	{
		count_vals[kk]=0;
		for (int ii=0; ii<MATSIZEX; ++ii)
		{
			for(int jj=0; jj< MATSIZEY; ++jj)
			{
				if(matrix[ii][jj] == kk)
					count_vals[kk]++;
			}
		}

	}
	}
}



void print_matrix(unsigned char matrix[MATSIZEX][MATSIZEY])
{
	for (int ii=0; ii<MATSIZEX; ++ii)
	{
		for(int jj=0; jj< MATSIZEY; ++jj)
		{
			std::cout << (int)matrix[ii][jj] << " ";
		}
		std::cout << std::endl;
	}

}

int main()
{
	
	unsigned char matrix[MATSIZEX][MATSIZEY];
	size_t count_vals[RANGE_HIGH];
	populate_matrix(matrix);
	count(matrix, count_vals);
	int sum=0;
	for(size_t index; index<RANGE_HIGH+1; index++)
	{
		//std::cout << index<< ": " << count_vals[index] << " " << std::endl;
		sum=sum+count_vals[index];
	}
	//std::cout << "SUM IS " << sum << " But should be "<< 128*128<<std::endl;

	
}


