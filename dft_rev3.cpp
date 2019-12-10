#include <iostream>
#include "complex.h"
#include "input_image.h"
#include <cmath>
#include <thread>

int num_threads = 8;

// Function to test thread call capability

void thread_call(int tid) 
{
     // std::cout << "Launched by thread " << tid << std::endl;
}


Complex DFT(int row,int k, int N, Complex* mat)
{
	float sumReal = 0.0f;
	float sumImag = 0.0f;
	float angle = 0.0f;
	Complex x = 0.0f;
	Complex newValue;
	Complex currValue;

	for(int n = 0; n < N; n++) // for every column in the row
	{
		angle = (2*M_PI*n*k) / N;
		//x = mat[row*N+n].mag(); //coefficient at the particular index in the matrix

		//sumReal += x.real*cos(angle);
		//sumImag -= x.real*sin(angle);

		/* new calc */
		currValue = mat[row*N+n];
		sumReal += currValue.real*cos(angle)+currValue.imag*sin(angle);
		sumImag += -1*currValue.real*sin(angle)+currValue.imag*cos(angle);

	}
	
	

	newValue.real = sumReal;
	newValue.imag = sumImag;

	return newValue;
}


void Trans(Complex* Mat, int N)
{
	Complex CopyMat[N*N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			CopyMat[i*N+j] = Mat[j*N+i]; // swap values
		}
	}

	for (int s = 0; s < N*N; s++)
	{
		Mat[s] = CopyMat[s]; // copy values back into the original Matrix
	}

	
}

int main(int argc, char* argv[])
{
	// retrieve input image, create image object
	InputImage image(argv[1]);
	int N = image.get_width(); // retrieve width from image

	// input values into matrix
	Complex* Matrix = image.get_image_data();
	Complex MatrixCopy[N*N];

	// define threads
	std::thread t[num_threads];

	// logic to determine how many threads to run //
	// if there are more rows than threads
	int rows_per_thread = N/num_threads;
	if ((rows_per_thread) > 0)
	{
		if (rows_per_thread == 1) // if we have an equal number of rows and threads
		{
			//Launch a group of threads
	 		for (int i = 0; i < num_threads; ++i) 
		     {
		         t[i] = std::thread(thread_call, i);
		     }	
		}
		// if even
		else if (rows_per_thread % 2 == 0)
		{
			//Launch a group of threads
     		for (int i = 0; i < num_threads; ++i) 
		     {
		         t[i] = std::thread(thread_call, i);
		     }
		}
		// if odd
		else 
		{
			//Launch a group of threads
			num_threads++;

 			for (int i = 0; i < num_threads+1; ++i) 
	     	{
	        	t[i] = std::thread(thread_call, i);
	     	}

		}
	}
	else // if there are more threads than rows
	{
		// divide threads by rows
		int new_num_threads = num_threads/N;
		// if even
		if (new_num_threads % 2 == 0)
		{
			//Launch a group of threads
 			for (int i = 0; i < new_num_threads; ++i) 
	     	{
	        	t[i] = std::thread(thread_call, i);
	     	}

	     	// update num_threads
	     	num_threads = new_num_threads;

		}
		// if odd
		else
		{
			//Launch a group of threads
 			for (int i = 0; i < new_num_threads+1; ++i) 
	     	{
	        	t[i] = std::thread(thread_call, i);
	     	}

	     	// update num_threads
	     	num_threads = new_num_threads+1;

		}
	}


	// row calculation
	for (int r = 0; r < N; r++)
	{
		// compute the DFT for every frequency
		for (int c = 0; c < N; c++)
		{
			MatrixCopy[r*N+c] = DFT(r,c,N,&Matrix[0]);
		}

	}

	/*
	// print out the DFT array
	for (int l = 0; l < N; l++)
	{
		for (int m = 0; m < N; m++)
		{
			std::cout << MatrixCopy[l*N+m] << "\t"; 
		}
		std::cout << std::endl;
	}
	*/
	

	Trans(&MatrixCopy[0],N);

	// column calculation
	
	for (int r = 0; r < N; r++)
	{
		// compute the DFT for every frequency
		for (int c = 0; c < N; c++)
		{
			Matrix[r*N+c] = DFT(r,c,N,&MatrixCopy[0]);
		}

	}
	

	/*
	for (int c = 0; c < N; c++)
	{
		for (int r = 0; r <= N*N-N; r+=N)
		{
			Matrix[r+c] = DFT_Col(r,c,N,&MatrixCopy[0]);
			//std::cout << "index beyotch: " << r+c << std::endl;
		}
	}
	*/

	Trans(&Matrix[0],N);


	
	
	

	//std::cout << "This is the version you are looking for!" << std::endl;

	//Join the threads with the main thread
     for (int i = 0; i < num_threads; ++i) {
         t[i].join();
     }

     // print out the DFT array
	for (int l = 0; l < N; l++)
	{
		for (int m = 0; m < N; m++)
		{
			std::cout << Matrix[l*N+m] << "\t"; 
		}
		std::cout << std::endl;
	}


	// end program
	return 0;
}