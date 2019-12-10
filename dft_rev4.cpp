#include <iostream>
#include "complex.h"
#include "input_image.h"
#include <cmath>
#include <thread>

int num_threads = 8;

// Function to test thread call capability
void thread_call(int tid) 
{
      std::cout << "Launched by thread " << tid << std::endl;
}

// function to print the current matrix
void Matrix_Print(Complex* mat, int N)
{
	
	for (int l = 0; l < N; l++)
		{
			for (int m = 0; m < N; m++)
			{
				std::cout << mat[l*N+m] << "\t\t"; 
			}
			std::cout << std::endl;
		}
	
		
}


// DFT calculation!
void DFT(int t_id, int t_length, int N, Complex* mat, Complex* matCopy)
{
	float sumReal = 0.0f;
	float sumImag = 0.0f;
	float angle = 0.0f;
	Complex x = 0.0f;
	Complex newValue;
	Complex currValue;
	int k = 0;
	int row = 0;


	int startPos = t_id*t_length;
	int endPos = startPos+t_length;

	//std::cout << "Thread No: " << t_id << "Index: " << startPos << std::endl;

	for (int i = startPos; i < endPos; i++) // iterate through the entire length of the thread
	{
		// reset values after each row calculation
	 	sumReal = 0.0f;
		sumImag = 0.0f;

		for(int n = 0; n < N; n++) // for every column in the row
		{
			k = i % N; // frequency to calculate is equal to the current index modulus row_width
			
			row = i/N; // find which row i'm on 

			angle = (2*M_PI*n*k) / N;
			
			/* new calc */
			currValue = mat[n+(row*N)]; // retreive the correct value at the right index from the main array

			sumReal += currValue.real*cos(angle)+currValue.imag*sin(angle); // sum the real portion
			sumImag += -1*currValue.real*sin(angle)+currValue.imag*cos(angle); // sum the imaginary portion

		}

		// assign values back into the matrix
		newValue.real = sumReal;
		newValue.imag = sumImag;

		matCopy[i] = newValue;
	}
	
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
	InputImage image(argv[2]);
	int N = image.get_width(); // retrieve width from image

	// input values into matrix
	Complex* Matrix = image.get_image_data();
	Complex MatrixCopy[N*N];

	// define threads
	std::thread t[num_threads];
	std::thread col_threads[num_threads];

	// logic to determine how many threads to run //
	// if there are more rows than threads
	int rows_per_thread = N/num_threads;
	if ((rows_per_thread) > 0)
	{
		if (rows_per_thread == 1) // if we have an equal number of rows and threads
		{
			// make no change to the amount of threads
			num_threads = num_threads;
		}
		// if even
		else if (rows_per_thread % 2 == 0)
		{
			// no change if there are an even number of rows
			num_threads = num_threads;
		}
		// if odd
		else 
		{
			// add one more thread if there are not an even number of rows for threads
			++num_threads;

		}
	}
	else // if there are more threads than rows
	{
		/*
		// divide threads by rows
		int new_num_threads = num_threads % N; // better algorithm for determining how many threads to spawn
		if (new_num_threads == 0)
		{
			new_num_threads = 1;
			num_threads = new_num_threads;
			rows_per_thread = N; // one thread needs to handle the whole matrix
		}
		// if even
		else if (new_num_threads % 2 == 0)
		{
	     	// update num_threads
	     	num_threads = new_num_threads;
	     	rows_per_thread = 

		}
		// if odd
		else
		{
	     	// update num_threads
	     	num_threads = new_num_threads+1;

		}
		*/
		num_threads = 1;
		rows_per_thread = N;
	}

	// define the block size for each thread
	int thread_length = rows_per_thread * N;


	// rows \\
	// thread - function call
	for (int i = 0; i < num_threads; ++i)
	{
		t[i] = std::thread(DFT,i,thread_length,N,&Matrix[0],&MatrixCopy[0]);
	}

  	//Join the threads with the main thread
     for (int i = 0; i < num_threads; ++i) 
     {
         t[i].join();
     }
	
	// transpose the matrix for simplicity
	Trans(&MatrixCopy[0],N);


	// columns \\
	// thread - function call
	for (int i = 0; i < num_threads; ++i)
	{
		col_threads[i] = std::thread(DFT,i,thread_length,N,&MatrixCopy[0],&Matrix[0]);
	}

	//Join the threads with the main thread
     for (int i = 0; i < num_threads; ++i) 
     {
         col_threads[i].join();
     }

	// transpose again for correct orientation
	Trans(&Matrix[0],N);

	// print out the DFT array
	Matrix_Print(&Matrix[0],N);

	
	// save to an output file
	//InputImage imageOut(argv[3]);
	image.save_image_data(argv[3],&Matrix[0],N,N);

	

	// end program
	return 0;
}