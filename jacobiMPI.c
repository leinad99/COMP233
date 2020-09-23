/************************************************
**  Daniel DeGraaf ~ COMP233 Fall 2019 ~ Capstone Project
**  Simulate the transfer of heat through a plate using OpenMP
**  Initial code taken from a government website and can be found at
**  (https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html)
**  compile with 'mpicc -o jacobiMPI.exe jacobiMPI.c -lm -openmp'
**  run with 'mpirun -machinefile ~/machines-openmpi -np 4 --bynode /tmp/node000-bccd/jacobiMPI.exe'
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include <omp.h>

const int WIDTH = 1000; //the width picture we're making 
const int HEIGHT = 1000; //height of the picture we're making
const float EPSILON = .01; //= .01, the lowest value of diffNorm before we break out of the loop
const int MAX_ITER = 500000; //Max number of times you can go through the loop

//temperature values for NSEW
const int NORTH = 100; 
const int SOUTH = 100;
const int EAST = 0;
const int WEST = 0;
const int MAXTEMP = 100; //max temperature that any node can be

const int MASTER = 0;

//Hybrid solution is *slightly* faster with 2 thread, this makes sense with the other resukts I have found in my testing
const int numThreads = 2; 

int main( argc, argv )
int argc;
char **argv;
{
    //Variable Dictionary
	int i; //used to initialize 2d array; 
	int r; //used to traverse rows 
	int c; //used to traverse columns

    int blueValue;
    int greenValue = 0; //only transition from red to blue
    int redValue;

	int iterationCount; //count how many iterations it takes to 
    float diffNorm; //use to know then the iteration has converged
	
	float **oldArray; //temperature values for each unit in the 2d array
	float **newArray; //used to store the new values of each unit from averaging the 4 adjecent indeces from oldArray

	FILE* fp; //pointer to file we will be writing to
	
	double start, end; //used to time how long the program takes
	int totalTime; //in seconds

    float globalDiffNorm; //global (accross all processors) value for the diffnorm
    MPI_Status status;
    int rank, size; //for MPI world
    int startIndex, endIndex; //start and and spots for each MPI processor

	int recvRank; //rank that the master rceives from to put the whole image together
	//1.0 Init -----------------------------------------------
    //MPI init 
    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

     if(rank == MASTER){
    	printf("\nDaniel DeGraaf ~ COMP233 Fall 2019 ~ Capstone Project \n\n");
        //start time
	    start = omp_get_wtime();
    }

    //Note that top and bottom processes have one less row of interior points
    startIndex = ((HEIGHT * rank) / size) ;
    endIndex = ((HEIGHT * (rank + 1)) / size);
    //so you don't touch the first row or the last row
    if (rank == 0){
        startIndex++; //so you don't touch the first row or the last row
    }        
    if (rank == size - 1){
        endIndex--;
    }

    /*
    **  I made a 1000 x 1000 array for all 4 processors instead of 250 x 1000 array
    **  for each processor and merge to one 1000 x 1000 array because from the testing 
    **  I and others did, doing it this way was ~24 seconds faster over the other method
    **  where each processor only has a 250 x 1000 array. I don't know why, but those were 
    **  my results and since I'm optimizing for speed this is what I choose to do. 
	*/
    //dynamically allocate the two 2D arrays
    oldArray = (float**) malloc(HEIGHT * sizeof(float*));
	for(i = 0; i < HEIGHT; i++) {
		oldArray[i] = (float*) malloc(WIDTH * sizeof(float));
	}
    newArray = (float**) malloc(HEIGHT * sizeof(float*));
	for(i = 0; i < HEIGHT; i++) {
		newArray[i] = (float*) malloc(WIDTH * sizeof(float));
	}

    //2.0 Process: Fill the arrays and loop through constantly changing the temperatures 
	//to an average of the indecies around it ------------------------------------------ 

    // Fill the data as specified
    for (r = 0; r < HEIGHT; r++) {
		for (c = 0; c < WIDTH; c++) {
			if(r == 0){ //NORTH
				oldArray[r][c] = NORTH;
				newArray[r][c] = NORTH;
			} else if(r == (HEIGHT - 1)) { //SOUTH
				oldArray[r][c] = SOUTH;
				newArray[r][c] = SOUTH;
			} else if(c == 0){//EAST
				oldArray[r][c] = EAST;
				newArray[r][c] = EAST;
			} else if(c == (WIDTH - 1)){ //WEST
				oldArray[r][c] = WEST;
				newArray[r][c] = WEST;
			} else{ //MIDDLE
				oldArray[r][c] = (NORTH + SOUTH + EAST + WEST) / 4.0;
				newArray[r][c] = (NORTH + SOUTH + EAST + WEST) / 4.0;
			}
		}
	}
	
    iterationCount = 0;
    do {
        //Share values accross MPI world 
        //Send the top up and the buttom down
       if (rank == MASTER) {
            //send the top up, and receive the button
            MPI_Send( oldArray[endIndex - 1], HEIGHT, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD );
            //receive the top from the one above you
            MPI_Recv( oldArray[endIndex], HEIGHT, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &status );
        } 
        else if (rank == size - 1){
            //receive your buttom from the one below you, then send your bottom down
            MPI_Recv( oldArray[startIndex - 1], HEIGHT, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status );
            MPI_Send( oldArray[startIndex], HEIGHT, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD );
        }  
        else{
            //send the top up
            MPI_Send( oldArray[endIndex - 1], HEIGHT, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD );
            MPI_Recv( oldArray[startIndex - 1], HEIGHT, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &status );
            //send the button down
            MPI_Send( oldArray[startIndex], HEIGHT, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD );
            MPI_Recv( oldArray[endIndex], HEIGHT, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &status );
        }

		//Compute new values but not on the boundary
		iterationCount++;
		diffNorm = 0.0f; //reset to 0

#pragma omp parallel for reduction(+:diffNorm) num_threads(numThreads) private(r, c) shared(newArray, oldArray)
		for (r = startIndex; r < endIndex; r++){
			for (c = 1; c < WIDTH - 1; c++) {
                //compute the new value from the surrounding nodes
				newArray[r][c] = (oldArray[r][c+1] + oldArray[r][c-1] + oldArray[r+1][c] + oldArray[r-1][c]) / 4.0f;
				diffNorm += (newArray[r][c] - oldArray[r][c]) * (newArray[r][c] - oldArray[r][c]);
			}
		}

		//put the new values into the old and make the pointer for the new point to the old 
		float **temp = oldArray;
        oldArray = newArray;
        newArray = temp;
        //merge the diffNorm values from all the processors
        MPI_Allreduce( &diffNorm, &globalDiffNorm, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD );

        //get the value for diffNorm
		globalDiffNorm = sqrtf(globalDiffNorm);

		//3.0 Output ---------------------------------------------
		if (iterationCount % 1000 == 0 && rank == MASTER){ //tell the user something is happening
			printf("At iteration %d, globalDiffNorm is %f\n", iterationCount, globalDiffNorm );
		} 
    } while (globalDiffNorm > EPSILON && iterationCount < MAX_ITER);

    //merge all the files together
    if (rank == MASTER) {
        //receive all the rows
        for(i = 1; i < size; i++){
            //+ WIDTH so we don't have blank line in the merged file
            MPI_Recv(oldArray[((HEIGHT * i) / size)], (WIDTH*HEIGHT)/size + WIDTH, MPI_FLOAT, i, 2, MPI_COMM_WORLD, &status ); 
        }
    }
    else{
        MPI_Send(oldArray[((HEIGHT * rank) / size)], (WIDTH*HEIGHT)/size + WIDTH, MPI_FLOAT, MASTER, 2, MPI_COMM_WORLD); 

    }
    //Put all the values in a ppm file
    if(rank == MASTER){
        //header for file
        fp = fopen("jacobiImage.ppm", "w+");
        fprintf(fp,"P3\n%d %d #image width (cols) and height (rows)\n", WIDTH, HEIGHT); //put the headers in the first line of the ppm 
        fprintf(fp, "#Daniel DeGraaf ~ COMP 233 ~ Laplace Heat Distribution\n");
        fprintf(fp, "#This image took %d iterations to converge\n", iterationCount);
        fprintf(fp, "255");
        for(r = 0; r < HEIGHT; r++) {
            for(c = 0; c < WIDTH; c++){
                if (c % 5 == 0){ //print 5 RGB values per line
                    fprintf(fp, "\n");
                }
                //use oldArray because it was last copied onto from newArray
                redValue = 255 * (oldArray[r][c] / MAXTEMP); //pixels at 100 degrees will be pure RED
                blueValue = 255 - (255 * (oldArray[r][c] / MAXTEMP)); //pixels at 0 degrees will be pure BLUE
                fprintf(fp, "%d %d %d ", redValue, greenValue, blueValue);
            }
        }
        //close file
        fclose(fp);
    }
    
    //4.0 Finish up ------------------------------------------
    //free up the data
    for(i = 0; i < WIDTH; i++){
        free(newArray[i]);
        free(oldArray[i]);
    }
	free(newArray);
    free(oldArray);
    
    if(rank == MASTER){
        //calculate the totalTime
        end = omp_get_wtime();
        totalTime = (int)(end - start); //total time is in seconds
        printf("\nFinal iteration count is %d\n", iterationCount);
        printf("Total runtime for the Hybrid solution [4 processors and 2 threads] is %d seconds.\n", totalTime);
        
        printf("\n\n\t<<< Normal Termination >>>\n\n");
    }
	
    MPI_Finalize();
    return 0;
}

