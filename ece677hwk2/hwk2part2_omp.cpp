#include <iostream>
#include <random>
#include <memory.h>
#include <omp.h>
#include <sys/time.h>

#define MATSIZEX 256
#define MATSIZEY 256
#define RANGE_HIGH 31
#define RANGE_LOW 0
#define WINDOWX 4
#define WINDOWY 4


/***********************************************
 * Populate a matrix of with a even distribution
 * of integers from RANGE_LOW to RANGE_HIGH
 **********************************************/
void populate_matrix(unsigned char matrix[MATSIZEX][MATSIZEY])
{
    std::random_device rd;
    std::mt19937 gen(rd());//mesenne twister engine
    std::uniform_int_distribution<unsigned char> dis(RANGE_LOW, RANGE_HIGH);
	
	//#pragma omp parallel
	//{
	//#pragma omp for
    for (int ii=0; ii<MATSIZEX; ++ii)
	{
		for(int jj=0; jj< MATSIZEY; ++jj)
		{
			matrix[ii][jj] = dis(gen);
		}
	}
	//}
}


long long getMicrotime()
{
        struct timeval currentTime;
        gettimeofday(&currentTime, NULL);
        return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}

void printtime(const char desc[10], int rank, int size, long long timer) 
{ 
        printf("omp> %s %i %i %f\n", desc, rank, size, (getMicrotime()-timer)/1e6); 
}

void print_matrix(unsigned char matrix[][256], int width, int height)
{
	for (int ii=0; ii<width; ++ii)
	{
		for(int jj=0; jj< height; ++jj)
		{
			printf("% 3d, ", matrix[ii][jj]);
		}
		printf("\n");
	}
	

}

/*************************************************
 *
 * Description: 
 * 	Take a x+MATSIZX by y+MATSIZEY slice
 * 	of a matrix and perform an RMS filter 
 * 	on that window. 
 *
 * **********************************************/
int filter_window(unsigned char matrix[MATSIZEX][MATSIZEY], size_t x, size_t y)
{
	double sum=0;

	for(int ii=x; ii<x+WINDOWX; ii++)
	{
		for(int jj=y; jj<y+WINDOWY; jj++)
		{
			sum = sum + matrix[ii][jj]*matrix[ii][jj];
		}
	}
	sum=sqrt(sum/(WINDOWX*WINDOWY));
	return (int) sum;
}

int main()
{
	
	unsigned char start_matrix[MATSIZEX][MATSIZEY];
	unsigned char final_matrix[MATSIZEX][MATSIZEY];

	//populate matrix randomly
	populate_matrix(start_matrix);

	//Deep copy the matrix
	//memcpy(final_matrix, start_matrix, MATSIZEX*MATSIZEY*sizeof(unsigned char) );
	int filter_val;

	long long timer = getMicrotime();
	#pragma omp parallel
	{
	printtime("init", omp_get_thread_num(), omp_get_num_threads(), timer );
	#pragma omp for
	for (int ii=0; ii<MATSIZEX-1; ++ii)
	{
		for(int jj=0; jj< MATSIZEY-1; ++jj)
		{
			filter_val = filter_window(start_matrix, ii, jj);
			for(int kk=1; kk<WINDOWX-1; kk++)
			{
				for(int ll=1; ll<WINDOWY-1; ll++)
				{
					final_matrix[ii+kk][jj+ll]=filter_val;
				}
			}
		
			//printf("_________________________________________________\n");
			//print_matrix(final_matrix, 8, 8);
			
		}
	}
	printtime("proc", omp_get_thread_num(), omp_get_num_threads(), timer);
	}
	
}


