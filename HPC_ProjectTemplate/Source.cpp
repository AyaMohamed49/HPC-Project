#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#pragma once

#include<mpi.h>
#include<stdio.h>

#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int width, height;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int *Red = new int[BM.Height * BM.Width];
	int *Green = new int[BM.Height * BM.Width];
	int *Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height*BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i*BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	width = BM.Width;
	height = BM.Height;
	return input;
}

void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i*width + j] < 0)
			{
				image[i*width + j] = 0;
			}
			if (image[i*width + j] > 255)
			{
				image[i*width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	//MyNewImage.Save("..//Data//Output//outputRes" + index + ".png");
	MyNewImage.Save("C:\\Users\\ayamo\\Documents\\GitHub\\HPC-Project\\Data\\OutPut\\10N\\" + index + ".png");
	cout << "result Image Saved " << index << endl;
}

int main()
{
	MPI_Init(NULL, NULL);

	int ImageWidth = 4, ImageHeight = 4, sizeOfImage = 256;
	int start_s, stop_s, TotalTime = 0;

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int s = 1; s <= 10; s++)
	{
		System::String^ imagePath;
		std::string img;
		img = "C:\\Users\\ayamo\\Documents\\GitHub\\HPC-Project\\Data\\Input\\10N\\"+to_string(s)+".png";

		imagePath = marshal_as<System::String^>(img);
		int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);

		int* pixel_intenisties = new int[256]{ 0 };
		double* probability = new double[256];
		double* comProbability = new double[256]{ 0 };
		int* floorComProbability = new int[256];
		/////////////////////////////////////////////////////////
		// parallel data
		int P_pixel_intenistiesLength = (width * height) / size - 1;

		//for step 1,2
		double* P_pixel_intenisties = new double[P_pixel_intenistiesLength];
		double* P_pixel_intenistiesArray = new double[256]{ 0 };
		int* sendCount1 = new int[size] {P_pixel_intenistiesLength};
		int* displs = new int[size] {0, P_pixel_intenistiesLength, P_pixel_intenistiesLength * 2, P_pixel_intenistiesLength * 3, P_pixel_intenistiesLength * 4, P_pixel_intenistiesLength * 5, P_pixel_intenistiesLength * 6, P_pixel_intenistiesLength * 7};
		int reminder = (width * height) % size - 1;
		sendCount1[size - 1] = reminder;

		//for step 3
		double* Final_P_pixel_intenistiesArray = new double[256]{ 0 };
		double localSum = 0;
		double* P_comutative_Probability = new double[256]{ 0 };

		//for step 4 and 5 
		double* Floor_P_comutative_Probability = new double[256]{ 0 };
		double* FinalOutputsmall = new double[P_pixel_intenistiesLength] {0};
		double* FinalOutput = new double[256]{ 0 };

		start_s = clock();

		/////////////////////////////////////////////////////////

		//sequential code
		if (size == 1)
		{
			if (rank == 0)
			{
				//step #1
				for (int i = 0; i < width * height; i++)
					pixel_intenisties[imageData[i]]++;

				//step #2
				for (int i = 0; i < sizeOfImage; i++)
					probability[i] = (double)pixel_intenisties[i] / (double)(width * height);

				//step #3
				double sum = 0;
				for (int i = 0; i < sizeOfImage; i++)
				{
					sum += probability[i];
					comProbability[i] = sum;
				}

				//step #4
				for (int i = 0; i < sizeOfImage; i++)
					floorComProbability[i] = floor(comProbability[i] * 256);

				//step #5
				for (int i = 0; i < width * height; i++)
					imageData[i] = floorComProbability[imageData[i]];
			}

		}
		//parallel code////////////////////////////
		else
		{


			//step 1
			cout << "FFF\n";
			MPI_Scatterv(imageData, sendCount1, displs, MPI_INT, P_pixel_intenisties, P_pixel_intenistiesLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//MPI_Scatter(imageData, P_pixel_intenistiesLength, MPI_DOUBLE, P_pixel_intenisties, P_pixel_intenistiesLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			for (int i = 0; i < P_pixel_intenistiesLength; i++)
			{
				P_pixel_intenistiesArray[(int)P_pixel_intenisties[i]]++;
			}
			MPI_Reduce(P_pixel_intenistiesArray, Final_P_pixel_intenistiesArray, P_pixel_intenistiesLength, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (rank == 0)
			{
				for (int i = 0; i < 256; i++)
				{
					cout << imageData[i] << " ";
				}
			}

			//step 3
			if (rank == 0)
			{
				for (int i = 0; i < 256; i++)
				{
					Final_P_pixel_intenistiesArray[i] /= (double)(width * height);
				}
				for (int i = 0; i < 256; i++)
				{
					localSum += Final_P_pixel_intenistiesArray[i];
					P_comutative_Probability[i] += localSum;
				}
				for (int i = 0; i < 256; i++)
				{
					Floor_P_comutative_Probability[i] = floor(P_comutative_Probability[i] * 255);
				}
			}

			//step 5 
			MPI_Scatterv(imageData, sendCount1, displs, MPI_INT, FinalOutputsmall, P_pixel_intenistiesLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//MPI_Scatter(imageData, P_pixel_intenistiesLength, MPI_DOUBLE, FinalOutputsmall, P_pixel_intenistiesLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			for (int i = 0; i < P_pixel_intenistiesLength; i++)
			{
				FinalOutputsmall[i] = Floor_P_comutative_Probability[(int)FinalOutputsmall[i]];
			}
			//MPI_Gather(FinalOutputsmall, P_pixel_intenistiesLength, MPI_DOUBLE, FinalOutput, P_pixel_intenistiesLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gatherv(FinalOutputsmall, P_pixel_intenistiesLength, MPI_DOUBLE, FinalOutput, sendCount1, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

			if (rank == 0)
			{
				createImage(imageData, ImageWidth, ImageHeight, 0);
			}
		}//HPC_ProjectTemplate

		
		////////////////////////////////////////////////////////

		
		start_s = clock();
		createImage(imageData, ImageWidth, ImageHeight, s);
		stop_s = clock();
		TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		cout << "time: " << TotalTime << endl;
		free(imageData);

	}
	MPI_Finalize();
	start_s = clock();
	stop_s = clock();
	TotalTime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
	cout << "Total time: " << TotalTime << endl;

	return 0;

}



