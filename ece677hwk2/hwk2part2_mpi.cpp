#include <iostream>
#include <random>
#include <memory.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>

#define THREAD_FINISHED -2
#define THREAD_SKIP -1
#define MATSIZEX 256
#define MATSIZEY 256
#define RANGE_HIGH 31
#define RANGE_LOW 0
#define WINDOWX 3
#define WINDOWY 3


/***********************************************
 * Populate a matrix of with a even distribution
 * of integers from RANGE_LOW to RANGE_HIGH
 **********************************************/
void populate_matrix(unsigned char **matrix)
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

long getMicrotime()
{
        struct timeval currentTime;
        gettimeofday(&currentTime, NULL);
        return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}



void print_matrix(unsigned char **matrix, int width, int height)
{
	for (int ii=0; ii<width; ++ii)
	{
		for(int jj=0; jj< height; ++jj)
		{
			if(jj!=height-1)
				printf("% 3d, ", matrix[ii][jj]);
			else
				printf("% 3d ", matrix[ii][jj]);
		}
		printf("\n");
	}
	

}

void print_matrix(const char *fname, unsigned char **matrix, int width, int height)
{
	FILE *fd=fopen(fname, "w");
	for (int ii=0; ii<width; ++ii)
	{
		for(int jj=0; jj< height; ++jj)
		{
			if(jj!=height-1)
				fprintf(fd,"% 3d, ", matrix[ii][jj]);
			else
				fprintf(fd,"% 3d ", matrix[ii][jj]);
		}
		fprintf(fd, "\n");
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

unsigned char **alloc_2d_char(int rows, int cols) 
{
    unsigned char *data = (unsigned char *)malloc(rows*cols*sizeof(unsigned char));
    unsigned char **array= (unsigned char **)malloc(rows*sizeof(unsigned char*));
	for (int i=0; i<rows; i++)
	{
		array[i] = &(data[cols*i]);
	}
    return array;
}


void printtime(const char desc[10], int rank, int size, long long timer) 
{ 
        printf("mpi> %s %i %i %f\n", desc, rank, size, (getMicrotime()-timer)/1e6); 
}




int main(int argc, char *argv[])
{
	long long timer;
	int rank, size;

	timer = getMicrotime();
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	char gomsg='G';
	char recvmsg;

	//unsigned char start_matrix[MATSIZEX][MATSIZEY]; 
	
	unsigned char ** start_matrix;
	unsigned char ** omatrix = alloc_2d_char( MATSIZEX, MATSIZEY );
	if(rank == 0)
		start_matrix = alloc_2d_char( MATSIZEX, MATSIZEY );
	
	if(rank == 0)
	{
		populate_matrix( start_matrix );
		//print_matrix( start_matrix, 8, 8 );
		//printf("and the num is %i\n", start_matrix[4][1]);

		
	}


	//number of 3by3 blocks in matrix
	int nblocks = (int) MATSIZEX*MATSIZEY;

	//number of bocks to be summed by each process
	//Not including 0 b/c thread 0 does the cooordination
	int blocks_per_proc = (int) nblocks/(size)-1;

	
	
	unsigned char final_matrix[MATSIZEX][MATSIZEY];
	MPI_Status stat;
	if (rank == 0)
	{
		//populate matrix randomly
		//populate_matrix( start_matrix );
		//print_matrix( start_matrix, 8, 8 );
		//Deep copy the matrix
		//memcpy( final_matrix, start_matrix, MATSIZEX*MATSIZEY*sizeof(unsigned char) );
		print_matrix("mat.dat", start_matrix, 256, 256 );
		for(int threadnum=1; threadnum<size; threadnum++)
		{
			MPI_Send( &(start_matrix[0][0]), MATSIZEX*MATSIZEY, MPI_BYTE, threadnum, 0, MPI_COMM_WORLD );
		}
	}

	//All threads recieve omatrix from thread0
	if(rank != 0)
	{
		MPI_Recv( &(omatrix[0][0]), MATSIZEX*MATSIZEY, MPI_BYTE, 0, 0, MPI_COMM_WORLD, &stat );
	}
		
	printtime("init", rank, size, timer);

	int startx=WINDOWX*(rank-1);
	int starty=0;
	int block_row=0;
	int block_col=0;
	double sum;

	//value and index of filter.
	int summsg[3];
	int sumrcv[3];
	int ndone=0;

	//thread0 processes responses
	if(rank == 0)
	{
		print_matrix("mat0.dat", start_matrix, MATSIZEX, MATSIZEY);
		timer = getMicrotime();
		for(;;)
		{
			MPI_Recv(sumrcv, 3, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
			if (sumrcv[0] == THREAD_FINISHED)
				ndone++;
			else if(sumrcv[0] != THREAD_SKIP)
			{
				start_matrix[sumrcv[1]][sumrcv[2]] = sumrcv[0];
			}
			if(ndone == size-1)//all threads have finished
				break;
		}
		printtime("Thread0", rank, size, timer );
		print_matrix( "mat1.dat", start_matrix, MATSIZEX, MATSIZEY );
	}


	else
	{//other threads work on the filter
		timer=getMicrotime();
		while(starty < MATSIZEY)
		{
			sum=0;
			if(startx > MATSIZEX)
			{//advance the column by 1
				startx = startx%MATSIZEX;
				block_row++;
				starty+=1;
			}

			//ignore the edges
			if(startx+WINDOWX > MATSIZEX || starty+WINDOWY > MATSIZEY)
			{
				summsg[0] = THREAD_SKIP;
				summsg[1] = startx;
				summsg[2] = starty;
				// we are at the edge so send skip message.
				
				MPI_Send(summsg, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
			else
			{
				if(rank == 1)
				for(int ii=startx; ii<startx+WINDOWX; ii++)
				{
					for(int jj=starty; jj<starty+WINDOWY; jj++)
					{
						sum+=omatrix[ii][jj]*omatrix[ii][jj];
					}
				}
				sum = sqrt(sum/(WINDOWX*WINDOWY));
				summsg[0] = (int) sum;
				summsg[1] = startx+1;
				summsg[2] = starty+1;
				//report the filtered value to thread0
				MPI_Send(summsg, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
			
			}
			//set startx for the next block to manipulate. 
			startx += 1*(size-1);
		}
		summsg[0] = THREAD_FINISHED;
		MPI_Send(summsg, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
		printtime("filter", rank, size, timer );
	}
	MPI_Finalize();
}


