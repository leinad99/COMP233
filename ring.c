/************************************************
** Daniel DeGraaf ~ COMP233 Fall 2019 ~ Ring around the Rosie
** See how fast it takes for a message to pass through all 16 cores
*************************************************/


#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
//function to print to a csv file
void createCSV(int results[12][4]);

int main(argc, argv)
int argc;
char** argv; {
	//variable dictionary
	int world_rank; //who am I
	int world_size; //how big is the world
	double* msg; //message each processor is passing to the next one
	double startTime; //time we start the clock
	double finishTime; //time we end the clock
	double totalTime; //how long it took to pass the message all the way aound
	double maxTime; //max time it took to pass the message all the way around
	double minTime; //min time it took to pass the message all the way around
	double avgTime; //avg time it took to pass the message all the way around
	int i; //iterator used to populate the message with values
	int timesThrough; //counter to keep track of how many times throuugh the program you've gone
	int messageSize; //size of the message for the current test
	int resultsCounter = 0; //counter so that you store the results in the correct spot
	double results[12][4]; //store 12 iterations for different message sizes and 4 values per iteration
	int col; //keep track of where to print the data to the file
	int row;
	FILE* fp; //pointer to file we will be writing to
	const int MAXMESSAGE = 1000000; //maximum message size for testing
	const int MAXITR = 10; //maximum times we will run for each message size to collect the data
	//We are setting master to 3 and not 0 because the way our cluster is set up, 
	//3 is the first node on the master machine, not 0
	//therefore in order to print to a csv file that we can access from master
	//node 3 must do all the calculations and print to the csv file
	const int MASTER = 3; 

	//1.0 Init -----------------------------------------------
	MPI_Init(&argc, &argv); //start the MPI world
	//Initialize Variables
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); //how many cores you have
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); //what core you are
	MPI_Status status;
	if (world_rank == 0) {//MASTER
		printf("\n\nDaniel DeGraaf ~ COMP233 Fall 2019 ~ Ring around the Rosie\n");
	}
	//make a message of 1000000 as that is the maximum amount of doubles you will ever put into it
	msg = (double*)malloc(1000000 * sizeof(double));

	//message sizes of 1000, 50000, 100000, 200000, 300000, ..., 1000000
	for (messageSize = 1000; messageSize <= MAXMESSAGE; messageSize += 100000) {
		//change message size for irregular jumps in message size at the beginning
		if (messageSize == 101000) {
			messageSize = 50000;
		}
		else if (messageSize == 150000) {
			messageSize = 100000;
		}

		maxTime = 0; //initialize them to 0
		avgTime = 0;
		minTime = 9999999; //set it equal to a bigger time then the message will ever reach

		for (timesThrough = 0; timesThrough < MAXITR; timesThrough++) {
			//2.0 Process: have the message pass from node to node until it
			//goes all the way through and back to the master
			if (world_rank == MASTER) { //MASTER
				for (i = 0; i < messageSize; i++) { //populate the message
					//fill the message with (sizeOfMsg + timeTestNum).  
					//The fifth test of the 1000-element message will have an array of all “1005.0” in it.
					msg[i] = messageSize + timesThrough + 1;
				}
				startTime = MPI_Wtime();
				MPI_Send(msg, messageSize + 1, MPI_DOUBLE, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD);
				MPI_Recv(msg, messageSize + 1, MPI_DOUBLE, (world_rank - 1) % world_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//validate that the message sent is the same as the message received.  
				if (msg[0] != (double)(messageSize + timesThrough + 1)) {
					printf("ERROR!!! Test message of size %d failed to pass correctly", messageSize);
				}
				finishTime = MPI_Wtime();
			}
			else { //SLAVE
				MPI_Recv(msg, messageSize + 1, MPI_DOUBLE, (world_rank - 1) % world_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(msg, messageSize + 1, MPI_DOUBLE, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD);
			}

			if (world_rank == MASTER) { //find the max, min, and average times
				totalTime = finishTime - startTime;
				avgTime += totalTime;
				if (totalTime > maxTime) {
					maxTime = totalTime;
				}
				if (totalTime < minTime) {
					minTime = totalTime;
				}
			}
		}
		//3.1 Output to console ---------------------------------------------
		if (world_rank == MASTER) { //MASTER
			avgTime = avgTime / 10; //to get the avg time around the whole loop
			//printf("\nData from message size of %d\n", messageSize);
			//printf("Max Time: %.5e seconds\n", maxTime);
			//printf("Min Time: %.5e seconds\n", minTime);
			//printf("Average Time: %.5e seconds\n", avgTime);
			results[resultsCounter][0] = (double) messageSize * 8; //multiply by 8 because each double is 8 bits
			//divide by 16 to get the time of one send and receive, not the whole loop
			results[resultsCounter][1] = minTime / 16; 
			results[resultsCounter][2] = maxTime / 16;
			results[resultsCounter][3] = avgTime / 16;
			resultsCounter++;
		}
	}
	//3.2 Output to file results.csv ---------------------------------------------
	if (world_rank == MASTER) { //MASTER
		//write the results to the file
		fp = fopen("results.csv", "w+"); 
		col = row = 0;
		fprintf(fp,"Message Length, Min Time, Max Time, Avg Time"); //put the headers in the first line of the csv
		for (col = 0; col < 12; col++) {
			fprintf(fp, "\n");//new line after every row
			for (row = 0; row < 4; row++) {
				if (row == 0) {//don't print  comma before the first result
					fprintf(fp, "%lf", results[col][row]);
				}
				else {
					fprintf(fp, ", %lf", results[col][row]);
				}
			}
		}
		//4.0 Finish up ------------------------------------------
		fclose(fp);
		printf("\n\n\t<<<Normal Termination>>>\n\n");
	}
	free(msg); //free the message on all the hosts
	MPI_Finalize();
	return 0;
}
