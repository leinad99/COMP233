/************************************************
** Daniel DeGraaf ~ COMP233 Fall 2019 ~ Jacobi OpenMP
** Simulate the transfer of heat through a plate using OpenMP
** Initial code taken from a government website and can be found at
** (https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src/jacobi/C/main.html)
** compile with: gcc -fopenmp jacobiOpenMP.c -o jacobiOpenMP.exe -lm
** run with: ./jacobiOpenMP.exe
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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

int numThreads = 4;

int main() { 
    //Variable Dictionary
	int i; //used to initialize 2d array; 
	int r; //used to traverse rows 
	int c; //used to traverse columns

    int blueValue;
    int greenValue = 0; //only transition from red to blue
    int redValue;

	int iterationCount; //count how many iterations it takes to 
    double diffNorm; //use to know then the iteration has converged
	
	float **oldArray; //temperature values for each unit in the 2d array
	float **newArray; //used to store the new values of each unit from averaging the 4 adjecent indeces from oldArray
    
	FILE* fp; //pointer to file we will be writing to
	
	double start, end; //used to time how long the program takes
	int totalTime; //in seconds
	
	//1.0 Init -----------------------------------------------
	printf("\nDaniel DeGraaf ~ COMP233 Fall 2019 ~ Jacobi \n\n");
	//start time
	start = omp_get_wtime();
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
			} else{
				oldArray[r][c] = (NORTH + SOUTH + EAST + WEST) / 4.0;
				newArray[r][c] = (NORTH + SOUTH + EAST + WEST) / 4.0;
			}
		}
	}
	
    iterationCount = 0;
    do {
		//Compute new values
		iterationCount++;
		diffNorm = 0.0; //reset to 0
#pragma omp parallel for reduction(+:diffNorm) num_threads(numThreads) private(r, c) shared(newArray, oldArray)
		for (r = 1; r < HEIGHT - 1; r++){ //start at 1 and end at HEIGHT - 1 because the boundary values don't change
			for (c = 1; c < WIDTH - 1; c++) {
				newArray[r][c] = (oldArray[r][c+1] + oldArray[r][c-1] + oldArray[r+1][c] + oldArray[r-1][c]) / 4.0;
				diffNorm += (newArray[r][c] - oldArray[r][c]) * (newArray[r][c] - oldArray[r][c]);
			}
		}

		//get the value for diffNorm
		diffNorm = sqrt(diffNorm);

		//put the new values into the old and make the pointer for the new point to the old 
		float **temp = oldArray;
        oldArray = newArray;
        newArray = temp;
		
		//3.0 Output ---------------------------------------------
		if (iterationCount % 1000 == 0){ //tell the user something is happening
			printf("At iteration %d, diffNorm is %e\n", iterationCount, diffNorm );
		} 
    } while (diffNorm > EPSILON && iterationCount < MAX_ITER);
	
    //Put all the values in a ppm file
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
    
    //4.0 Finish up ------------------------------------------
    //close file
    fclose(fp);

    //free up the data
	for(i = 0; i < (HEIGHT); i++) {
		free(oldArray[i]);
	}
	free(oldArray);
	for(i = 0; i < (HEIGHT); i++) {
		free(newArray[i]);
	}
	free(newArray);

	//calculate the totalTime
	end = omp_get_wtime();
	totalTime = (int)(end - start); //total time is in seconds
    printf("\nFinal iteration count is %d, with the diffNorm being %e\n", iterationCount, diffNorm );
	printf("Total runtime for the OpenMP solution with %d threads is %d seconds.\n", numThreads, totalTime); //in minutes:seconds 
	
	printf("\n\n\t<<< Normal Termination >>>\n\n");
    return 0;
}

