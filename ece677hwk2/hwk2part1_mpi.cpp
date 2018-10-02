#include <iostream>
#include <random>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#define MATSIZEX 128
#define MATSIZEY 128
#define RANGE_HIGH 31
#define RANGE_LOW 0
#include <sys/time.h>

void populate_matrix(unsigned char matrix[MATSIZEX][MATSIZEY])
{
    std::random_device rd;
    std::mt19937 gen(rd());//mesenne twister engine
    std::uniform_int_distribution<unsigned char> dis(RANGE_LOW, RANGE_HIGH);
	

    for (int ii=0; ii<MATSIZEX; ++ii)
	{
		for(int jj=0; jj< MATSIZEY; ++jj)
		{
			matrix[ii][jj] = dis(gen);
		}
	}
}


/**
 * Returns the current time in microseconds.
 */
long getMicrotime()
{
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}

void inline count( unsigned char matrix[MATSIZEX][MATSIZEY], size_t count_vals[] )
{
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

void inline count_onenum( unsigned char matrix[MATSIZEX][MATSIZEY], size_t count_vals[], int num )
{
		count_vals[num]=0;
		for (int ii=0; ii<MATSIZEX; ++ii)
		{
			for(int jj=0; jj< MATSIZEY; ++jj)
			{
				if(matrix[ii][jj] == num)
					count_vals[num]++;
			}
		}

}

/***********************************************************
 * Count values in partitions of the matrix
*********************************************************/
void inline count_chunk( unsigned char matrix[MATSIZEX][MATSIZEY], int count_vals[], size_t startx, size_t starty, size_t chunksizex, size_t chunksizey)
{
	assert((startx+chunksizex) <= MATSIZEX );
	assert((starty+chunksizey) <= MATSIZEY );
	for(int kk=RANGE_LOW; kk<RANGE_HIGH+1; kk++)
	{
		count_vals[kk]=0;
		for (int ii=startx; ii<startx+chunksizex; ++ii)
		{
			for(int jj=starty; jj<starty+chunksizey; ++jj)
			{
				if(matrix[ii][jj] == kk)
				{
					count_vals[kk]++;
				}
			}
		}

	}
}

void print_matrix(unsigned char matrix[][MATSIZEY], int startx, int starty, int width, int height)
{
	for (int ii=startx; ii<startx+width; ++ii)
	{
		for(int jj=starty; jj< starty+height; ++jj)
		{
			printf("% 3d, ", matrix[ii][jj]);
		}
		printf("\n");
	}
	printf("\n");
	

}

void print_count(int count_vals[])
{
	int sum;
	for(size_t index=0; index<RANGE_HIGH+1; index++)
	{
		std::cout << index<< ": " << count_vals[index] << " ";
		sum=sum+count_vals[index];
	}
	std::cout << std::endl << "SUM IS " << sum << " But should be "<< MATSIZEX*MATSIZEY<<std::endl;
}

void printtime(const char desc[10], int rank, int size, long long timer)
{
	printf("mpi> %s %i %i %f\n", desc, rank, size, (getMicrotime()-timer)/1e6);
}

int main(int argc, char * argv[])
{
	unsigned char matrix[MATSIZEX][MATSIZEY];
	memset(matrix, 0, sizeof(char)*MATSIZEX*MATSIZEY);
	populate_matrix(matrix);
	exit(0);
	int rank, size, ii, jj;
	long long timer;
	long long timers[4];

	timer = getMicrotime();
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printtime("init", rank, size, timer);
	//printf("mpi> init %i %i %f\n", rank, size, (getMicrotime()-timer)/1e6);

	int local_count_vals[RANGE_HIGH+1];
	int global_count_vals[RANGE_HIGH+1];
	int fudge_factor=0;
	int startx, starty, width, height, sides;

	timer=getMicrotime();//restart timer
	for(int ii=0; ii<=RANGE_HIGH; ii++)
	{
		local_count_vals[ii]=0;
		global_count_vals[ii]=0;
	}

	
	

	//Break the matrix into partitions based on the rank of the thread. 
	count_chunk(matrix, local_count_vals, (rank)*MATSIZEX/size, 0, MATSIZEX/size, MATSIZEY );
	
	printtime("calc", rank, size, timer);
	/*printf("rank: %d, chunkx: %d-%d, chunky: %d-%d\n", 
			rank, 
			(rank)*MATSIZEX/size, 
			(rank)*MATSIZEX/size+MATSIZEX/size,  
			(rank)*MATSIZEY/size, 
			(rank)*MATSIZEY/size+MATSIZEX/size);
	*/

	timer=getMicrotime();//restart timer
	MPI_Reduce(local_count_vals, global_count_vals, RANGE_HIGH+1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	printtime("comm", rank, size, timer);

	int sum=0;

	if(rank == 0)
	{
		//print_matrix(matrix, 0,0,MATSIZEX, MATSIZEY);
		
		for(size_t index=0; index<RANGE_HIGH+1; index++)
		{
			//std::cout << index<< ": " << global_count_vals[index] << " ";
			sum=sum+global_count_vals[index];
		}
		//std::cout << std::endl << "SUM IS " << sum << " But should be "<< MATSIZEX*MATSIZEY<<std::endl;
		
		//print_count(global_count_vals);
	}

	timers[rank] = getMicrotime()-timers[rank];

	timer=getMicrotime();//restart timer
	MPI_Finalize();
	printtime("final", rank, size, timer);
}


