#include <iostream>
#include <random>
#include <omp.h>
#include <sys/time.h>

#define MATSIZEX 128
#define MATSIZEY 128
#define RANGE_HIGH 31
#define RANGE_LOW 0

void populate_matrix(unsigned char matrix[MATSIZEX][MATSIZEY])
{
    std::random_device rd;
    std::mt19937 gen(rd());//mesenne twister engine
    std::uniform_int_distribution<unsigned char> dis(RANGE_LOW, RANGE_HIGH);
	
	//omp_set_num_threads(2);
    for (int ii=0; ii<MATSIZEX; ++ii)
	{
		for(int jj=0; jj< MATSIZEY; ++jj)
		{
			matrix[ii][jj] = dis(gen);
		}
	}
}


long getMicrotime()
{
        struct timeval currentTime;
        gettimeofday(&currentTime, NULL);
        return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}

void printtime(const char desc[10], int rank, int size, long long timer) 
{ 
        printf("omp> %s %i %i %f\n", desc, rank, size, (getMicrotime()-timer)/1e6); 
}

void inline count( unsigned char matrix[MATSIZEX][MATSIZEY], size_t count_vals[] )
{
	long long timer = getMicrotime();
	#pragma omp parallel
	{
	printtime("init", omp_get_thread_num(), omp_get_num_threads(), timer );
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
	printtime("count", omp_get_thread_num(), omp_get_num_threads(), timer );
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

	for(size_t index=0; index<RANGE_HIGH+1; index++)
	{
		//std::cout << index<< ": " << count_vals[index] << " " << std::endl;
		//sum=sum+count_vals[index];
	}
	std::cout << "checksum: SUM IS " << sum << " But should be "<< MATSIZEX*MATSIZEY<<std::endl;
}


