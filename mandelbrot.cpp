/************************************************
** Daniel DeGraaf ~ COMP233 Fall 2019 ~ Multicolored Mandelbrot Image
** Creates a multicolored mandlebrot image with parallel computing
** Credit to Farouk Ounanes for the original code
** (https://github.com/ChinksofLight/mandelbrot_cpp/blob/master/mandelbrot.cpp)
*************************************************/
#include <iostream>
#include <fstream>
#include <complex>
#include <omp.h>

using namespace std;

const int WIDTH = 3000;
const int HEIGHT = 3000; //outputs a 3000px by 3000px image 
const int MAX_ITR = 300; //Given max number of iterations
int numThreads = 4; //adjust thread count if your machine has more cores

//Struct to hold the red, green, and blue values for a RBG pixel
struct RGB { //default color is black
	int red = 0;
	int green = 0;
	int blue = 0;
};

int value(int x, int y); 
void populateRGBValsArray(RGB RGBvalues[MAX_ITR], int values[WIDTH * HEIGHT]);

int main() {
	//Declare variables
	int* values = new int[WIDTH * HEIGHT]; //value for each pixel calculated by value()
	RGB RGBvals[MAX_ITR]; //RGB value for each possible value in values[]
	ofstream my_Image("mandelbrot.ppm");
	stringstream concatString; 
	//used to store the string that will be eventually inputted into the file
	//Did this to save on the extra time to input to a file, by storing the whole file as
	//a string and then inputing it all to the file once.
	

	//1.0 Init -----------------------------------------------
	cout << "Daniel DeGraaf ~ COMP233 Fall 2019 ~ Multicolored Mandelbrot Image\n";
	cout << "\tStores a multicolored mandelbrot image in a ppm format called mandelbrot.ppm \n\n";
	//Print 'UI' loading animation
	cout << "                                                 \\/ - Done\n";
	cout << "Loading";

	//2.0 Process: Calculate are the RGB values in a Parallelized for loop -------------
	populateRGBValsArray(RGBvals, values);
	
	//3.0 Output ---------------------------------------------
	if (my_Image.is_open()) {
		//Put in the needed info and heading
		concatString << "P3\n" << WIDTH << ' ' << HEIGHT << " 255\n";
		concatString << "#COMP 233 - Multicolored Mandelbrot Image\n";
		concatString << "#Author: Daniel DeGraaf, with starter code taken from Farouk Ouneanes\n";
		for (int i = 0; i < (WIDTH * HEIGHT); i++) {
			concatString << RGBvals[values[i]].red << ' ' << RGBvals[values[i]].green << ' ' << RGBvals[values[i]].blue;
			if (i % 5 == 0) //print 5 RGB values per line
				concatString << "\n";
			else
				concatString << ' ';
			if (i % 300000 == 0) //print loading dot
				cout << ".";
		}
		//convert stringstream to string and put that string in the file
		my_Image << concatString.str();
		
		//4.0 Finish up ------------------------------------------
		cout << "\nFinished!";
		//cleanup
		delete[] values;
		my_Image.close();

		cout << "\n\n\t<<< Normal Termination >>>\n";
		return 0;
	}
	else { //couldn't open the file
		cout << "\n\n\tError: Could not open the file" << endl;
		return -1;
	}
}

//Returns the escape iteration value (0-299)
int value(int x, int y) { 
	//Variable Dictionary
	complex<float> point((float)x / WIDTH - 1.5, (float)y / HEIGHT - 0.5);
	//we divide by the image dimensions to get values smaller than 1
	//then apply a translation
	complex<float> z(0, 0);
	int nb_iter = 0; //how many iterations it has gone through

	while (abs(z) < 2 && nb_iter < MAX_ITR) { //(MAX_ITR - 1) since you add one in the for loop 
		z = z * z + point;					  //and you don't want a value more than 299
		nb_iter++;
	}
	//Set the RGB value according to the number of iterations it passed through 
	if (nb_iter == MAX_ITR) 
		return 0; //0 represents the black inner circles of the image
	else 
		return nb_iter;
}

void populateRGBValsArray(RGB RGBvals[MAX_ITR], int values[WIDTH * HEIGHT]) {
	//Variable Dictionary
	int count[MAX_ITR] = { 0 };//set all the counts to 0
	int histogramVals[MAX_ITR] = { 0 }; //the histogram values for each of the counts, start them at 0
	double uDist[MAX_ITR]; //the u-Distrubution values for each of the histograms

	//populate the count and values array in parallel
	omp_set_num_threads(numThreads);
#pragma omp parallel
	{
		//Declare Variables
		int myID = omp_get_thread_num();
		int realNumThreads = omp_get_num_threads();

		int myStart = (long)(myID * HEIGHT) / realNumThreads;
		int myStop = (long)((myID + 1) * HEIGHT) / realNumThreads;

		for (int row = myStart; row < myStop; row++) {
			for (int col = 0; col < WIDTH; col++) {
				//Calculate the RGB value
				int val = value(row, col);
				//add one to how many times we have seen that iteration value appear 
				count[val] = count[val] + 1; 
				//set put the value in the right spot of values
				values[(row * WIDTH) + col] = val;
			}
			if(row % 300 == 0) 
				cout << ".";//print loading dot
		}
	}

	//set the histogram values
	for (int i = 1; i < MAX_ITR; i++) //don't calculate how many 0's there are
		//add the current count to the previous histogram
		histogramVals[i] = histogramVals[i - 1] + count[i]; 
	cout << ".";//print loading dot

	//populate the u-distribution
	for (int i = 0; i < MAX_ITR; i++) //don't calculate how many 0's there are
		uDist[i] = (float)histogramVals[i] / (float)histogramVals[MAX_ITR - 1];
	cout << ".";//print loading dot

	//populate the RGB vals for each uDist
	for (int i = 0; i < MAX_ITR; i++) {//don't calculate how many 0's there are
		RGB temp;
		//COLOR TEMPLATE
		//In this case do red, pink, and blue
		RGB white = { 255, 255, 255 };
		RGB black = { 0, 0, 0 };
		RGB color1 = { 191, 0, 64 };//Red: 191, 0, 64
		RGB color2 = { 128, 0, 128 };//Pink: 128, 0, 128
		RGB color3 = { 89, 0, 166 };//Blue: 89, 0, 166
		if (uDist[i] == 0)
			temp = black; 
			//this is to make the inner cicles black, we do this because in the value function
			//we set the val to 0 if the iterations equaled or exceeded the max
		else if (uDist[i] < .125) {//white -> black
			temp.red = white.red + (-white.red * (uDist[i] / .125));
			temp.green = white.green + (-white.green * (uDist[i] / .125));
			temp.blue = white.blue + (-white.blue * (uDist[i] / .125));
		}
		else  if (uDist[i] < .25) {//black -> color1
			temp.red = black.red + (color1.red * ((uDist[i] - .125) / .125));
			temp.green = black.green + (color1.green * ((uDist[i] - .125) / .125));
			temp.blue = black.blue + (color1.blue * ((uDist[i] - .125) / .125));
		}
		else  if (uDist[i] < .375) {//color1 -> black
			temp.red = color1.red + (-color1.red * ((uDist[i] - .25) / .125));
			temp.green = color1.green + (-color1.green * ((uDist[i] - .25) / .125));
			temp.blue = color1.blue + (-color1.blue * ((uDist[i] - .25) / .125));
		}
		else  if (uDist[i] < .5) {//black -> color2
			temp.red = black.red + (color2.red * ((uDist[i] - .375) / .125));
			temp.green = black.green + (color2.green * ((uDist[i] - .375) / .125));
			temp.blue = black.blue + (color2.blue * ((uDist[i] - .375) / .125));
		}
		else  if (uDist[i] < .625) {//color2 -> black
			temp.red = color2.red + (-color2.red * ((uDist[i] - .5) / .125));
			temp.green = color2.green + (-color2.green * ((uDist[i] - .5) / .125));
			temp.blue = color2.blue + (-color2.blue * ((uDist[i] - .5) / .125));
		}
		else  if (uDist[i] < .75) {//black -> color3
			temp.red = black.red + (color3.red * ((uDist[i] - .625) / .125));
			temp.green = black.green + (color3.green * ((uDist[i] - .625) / .125));
			temp.blue = black.blue + (color3.blue * ((uDist[i] - .625) / .125));
		}
		else if (uDist[i] < .875) {//color3 -> black
			temp.red = color3.red + (-color3.red * ((uDist[i] - .75) / .125));
			temp.green = color3.green + (-color3.green * ((uDist[i] - .75) / .125));
			temp.blue = color3.blue + (-color3.blue * ((uDist[i] - .75) / .125));
		}
		else if (uDist[i] < 1.0) {//black -> white
			temp.red = black.red + (white.red * ((uDist[i] - .875) / .125));
			temp.green = black.green + (white.green * ((uDist[i] - .875) / .125));
			temp.blue = black.blue + (white.blue * ((uDist[i] - .875) / .125));
		}
		RGBvals[i] = temp;
	}
	cout << ".";//print loading dot
}