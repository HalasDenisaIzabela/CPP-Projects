// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "common.h"
#include <algorithm>
#include <queue>
#include <random>
using namespace std;
using namespace cv;

void testOpenImage()
{
	char fname[MAX_PATH];
	while(openFileDlg(fname))
	{
		Mat src;
		src = imread(fname);
		imshow("image",src);
		waitKey();
	}
}

void MyCallBackFunc(int event, int x, int y, int flags, void* param)
{
	//More examples: http://opencvexamples.blogspot.com/2014/01/detect-mouse-clicks-and-moves-on-image.html
	Mat* src = (Mat*)param;
	if (event == EVENT_LBUTTONDOWN)
		{
			printf("Pos(x,y): %d,%d  Color(RGB): %d,%d,%d\n",
				x, y,
				(int)(*src).at<Vec3b>(y, x)[2],
				(int)(*src).at<Vec3b>(y, x)[1],
				(int)(*src).at<Vec3b>(y, x)[0]);
		}
}

void showHistogram(const std::string& name, int* hist, const int  hist_cols, const int hist_height)
{
	Mat imgHist(hist_height, hist_cols, CV_8UC3, CV_RGB(255, 255, 255)); // constructs a white image

	//computes histogram maximum
	int max_hist = 0;
	for (int i = 0; i < hist_cols; i++)
		if (hist[i] > max_hist)
			max_hist = hist[i];
	double scale = 1.0;
	scale = (double)hist_height / max_hist;
	int baseline = hist_height - 1;

	for (int x = 0; x < hist_cols; x++) {
		Point p1 = Point(x, baseline);
		Point p2 = Point(x, baseline - cvRound(hist[x] * scale));
		line(imgHist, p1, p2, CV_RGB(255, 0, 255)); // histogram bins colored in magenta
	}

	imshow(name, imgHist);
}

void showHistogram(const std::string& name, double* hist, const int  hist_cols, const int hist_height)
{
	Mat imgHist(hist_height, hist_cols, CV_8UC3, CV_RGB(255, 255, 255)); // constructs a white image

	//computes histogram maximum
	double max_hist = 0;
	for (int i = 0; i < hist_cols; i++)
		if (hist[i] > max_hist)
			max_hist = hist[i];
	double scale = 1.0;
	scale = (double)hist_height / max_hist;
	int baseline = hist_height - 1;

	for (int x = 0; x < hist_cols; x++) {
		Point p1 = Point(x, baseline);
		Point p2 = Point(x, baseline - cvRound(hist[x] * scale));
		line(imgHist, p1, p2, CV_RGB(255, 0, 255)); // histogram bins colored in magenta
	}

	imshow(name, imgHist);
}



//////////// FUNCTII ////////////////
void testOpenImagesFld()
{
	char folderName[MAX_PATH];
	if (openFolderDlg(folderName) == 0)
		return;
	char fname[MAX_PATH];
	FileGetter fg(folderName, "bmp");
	while (fg.getNextAbsFile(fname))
	{
		Mat src;
		src = imread(fname);
		imshow(fg.getFoundFileName(), src);
		if (waitKey() == 27) //ESC pressed
			break;
	}
}

/*int histograme[][256] = { 0 };

void testOpenImagesFldModif()
{
	int index = 0;
	char folderName[MAX_PATH];
	if (openFolderDlg(folderName) == 0)
		return;
	char fname[MAX_PATH];
	FileGetter fg(folderName, "bmp");
	while (fg.getNextAbsFile(fname))
	{
		Mat src;
		src = imread(fname);
		int height = src.rows;
		int width = src.cols;

		//int hist[256] = { 0 };
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				uchar pixel = src.at<uchar>(i, j);
				//hist[pixel]++;
				histograme[index][pixel]++;
			}
		}
		cout << index << " ";
		imshow(fg.getFoundFileName(), src);
		showHistogram("histograma", histograme[index], 255, 255);
		index++;
		for (int i = 0; i < 256; i++) {
			histograme[index][i] = { 0 };
		}
		
		if (waitKey() == 27) //ESC pressed
			break;
	}
}*/

//////////////////// PROIECT //////////////////////////

//Aceasta functie face conversia in HSI pentru o imagine si returneaza media Hue
float convHSI(Mat src) {

	float media;
	Mat hsi = Mat(src.rows, src.cols, CV_8UC3);
	int r, g, b;
	double h, s, in;
	float sum = 0.;
	int index = 0;
	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			r = src.at<Vec3b>(i, j)[2];
			g = src.at<Vec3b>(i, j)[1];
			b = src.at<Vec3b>(i, j)[0];
				
			in = (b + g + r) * 1.00 / 3;

			int minimRGB = 0;
			minimRGB = min(r, min(b, g));
			if (in == 0)
				s = 0;
			else
				s = 1 - (minimRGB / in);

			if (s != 0)
			{
				h = 0.5 * ((r - g) + (r - b)) / sqrt(((r - g) * (r - g)) + ((r - b) * (g - b)));
				h = acos(h);

				if (b > g)
				{
					h = ((360 * PI) / 180.0) - h;
				}
				else 
					h = h;
			}
			else {
				h = 0;
			}
			hsi.at<Vec3b>(i, j)[0] = ((h * 180) / PI)*255/360;
			hsi.at<Vec3b>(i, j)[1] = s * 100;
			hsi.at<Vec3b>(i, j)[2] = in;
			if (hsi.at<Vec3b>(i, j)[1] != 0)
			{
				sum += h;
				index++;
			}
		}
	} 
	if (index != 0) {
		media = sum / (float)index;
	}
	else
		media = -1;
	return media;
}

/// Functia in care impart imaginea in bucati
void gridImagine(Mat src, Mat bucati[], int &index, int N, float V[], int &indexV)
{
	index = 0;
	Mat image_copy = src.clone();
	int height = src.rows;
	int width = src.cols;
	int nOK1 = height / N;
	int nOK2 = width / N;
	int x1 = 0;
	int y1 = 0;
	indexV = 0;
	for (int y = 0; y <= height; y = y + nOK1)
	{
		for (int x = 0; x <= width; x = x + nOK2)
		{
			Mat tiles;
			if ((height - y) < nOK1 || (width - x) < nOK2)
			{
				break;
			}
			y1 = y + nOK1;
			x1 = x + nOK2;
			if (x1 >= width && y1 >= height)
			{
					x1 = width - 1;
					y1 = height - 1;
					tiles = image_copy(Range(y, height), Range(x, width));
			}
			else if (y1 >= height)
			{
					y1 = height - 1;
					tiles = image_copy(Range(y, height), Range(x, x + nOK2));
			}
			else if (x1 >= width)
			{
					x1 = width - 1;
					tiles = image_copy(Range(y, y + nOK1), Range(x, width));
			}
			else
			{
					tiles = image_copy(Range(y, y + nOK1), Range(x, x + nOK2));
			}
			bucati[index] = tiles;
			index++;
			if (convHSI(tiles) != -1) {
				V[indexV] = convHSI(tiles);
				indexV++;
			}
		}
	}
}
/*void testGridImagine() {
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src = imread(fname, IMREAD_COLOR);
		Mat bucati[1000] = {};
		int index = 0;
		int N = 5;
		float V[1000] = {};
		int indexV = 0;
		gridImagine(src, bucati, index, N, V, indexV);
		imshow("Patched Image", src);
		waitKey();
	}
}*/

/// Functia in care pun in vector de Mat toate pozele dintr-un folder
void testOpenImagesFldModif2(Mat imagini[], int &index)
{
    index = 0;
	char folderName[MAX_PATH] = "C:\\Users\\User\\Desktop\\OpenCVApplication-VS2019_OCV451_basic\\OpenCVApplication-VS2019_OCV451_basic\\Images\\Proiect";
	if (openFolderDlg(folderName) == 0)
		return;
	char fname[MAX_PATH];
	FileGetter fg(folderName, "bmp");
	while (fg.getNextAbsFile(fname))
	{
		Mat src;
		src = imread(fname);
		imagini[index] = src;
		index++;
		if (waitKey() == 27) //ESC pressed
			break;
	}
}

/// Se ia fiecare poza din folder si se imparte in N*N bucati
/*void prelucrareImgFld() 
{
	Mat imagini[1000];
	int index = 0;
	int N;
	cin >> N;
	float V[100];
	int indexV;
	testOpenImagesFldModif2(imagini, index);
	for (int i = 0; i < index ; i++)
	{
		string num = "img";
		num.append(to_string(i));
		//imshow(num, imagini[i]);
		int index2 = 0;
		Mat bucati[1000] = {};
		gridImagine(imagini[i], bucati, index2, N, V, indexV);
		
		//aici s-a testat daca se face impartirea corect in bucati a fiecarei poze din folder ul ales
		for (int i = 0; i < index2; i++)
		{
			string num2 = "bucata";
			num2.append(to_string(i));
			num2.append("_");
			num2.append(num);
			cout << num2 << endl;
			imshow(num2, bucati[i]);
			
		}
		cout << endl;

		for (int i = 0; i < indexV; i++)
		{
			cout << V[i] << " ";
		}
		cout << num << endl;

	}
	waitKey();

}*/

/*void testGrid()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src = imread(fname, IMREAD_COLOR);
		Mat bucati[1000];
		int index = 0;
		int indexV;
		int N;
		cin >> N;
		float V[1000];
		gridImagine(src, bucati, index, N, V, indexV);
		for (int i = 0; i < index; i++)
		{
			string num = "bucata";
			num.append(to_string(i));
			//imshow(num, bucati[i]);
			//cout << num << endl;
		}
		imshow("Imagine", src);

    	waitKey();
	}
}*/

float similitudini(float Voriginal[], float V[], int N)
{
	float k = 0;
	int x = 0;
	for (int i = 0; i < N * N; i++)
	{
		if (Voriginal[i] >= 0 && V[i] >= 0)
		{
			float minimum = min(abs(Voriginal[i] - V[i]), 2 * PI - abs(Voriginal[i] - V[i]));
			k += minimum * minimum;
			x++;
		}
	}
	float C = PI * sqrt(x);
	float D = (1 / C) * sqrt(k);
	float S;
	if (x == 0)
	{
		 S = 0;
	}
	else if(x>0)
	{
		S = 1 - D;
	}
	return S;
}

/*void proiect(Mat similare[])
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src = imread(fname, IMREAD_COLOR);
		Mat imagini[1000];
		Mat bucatiOriginal[1000];
		int indexOriginal = 0;
		float Voriginal[1000];
		int indexVOriginal;
		int index = 0;
		int N;
		cin >> N;
		float V[1000];
		int indexV;
		float simili[1000] = {0.0};
		gridImagine(src, bucatiOriginal, indexOriginal, N, Voriginal, indexVOriginal);
		testOpenImagesFldModif2(imagini, index);
		for (int i = 0; i < index; i++)
		{
			//string num = "img";
			//num.append(to_string(i));
			//imshow(num, imagini[i]);
			int index2 = 0;
			Mat bucati[1000] = {};
			gridImagine(imagini[i], bucati, index2, N, V, indexV);
			for (int j = 0; j < N*N; j++)
			{
				simili[i] = similitudini(Voriginal, V, N);
			}
		}
		float similiCop[100];
		for (int i = 0; i < index; i++)
		{
			similiCop[i] = simili[i];
		}
		float temp = 0.0f;
		for (int i = 0; i < index; i++)
		{
			for (int j = 0; j < index - i - 1; j++)
			{
				if (similiCop[j] > similiCop[j + 1])
				{
					temp = similiCop[j + 1]; 
					similiCop[j + 1] = similiCop[j];
					similiCop[j] = temp;
				}
			}
		}
		for (int i = 0; i < index; i++)
		{
			//cout << "img" << i << "___" << simili[i]<<endl;
		}
		for (int i = 0; i < index; i++)
		{
			//cout << "imgcop" << i << "___" << similiCop[i]<<endl;
		}
		for (int i = 0; i < index; i++)
		{
			if (simili[i] == similiCop[index - 1])
			{
				//cout << endl;
				//cout << "img" << i << "___" << simili[i]<<endl;
				//imshow(to_string(i), imagini[i]);
				similare[0] = imagini[i];

			}
			if (simili[i] == similiCop[index - 2])
			{
				//cout << endl;
				//cout << "img" << i << "___" << simili[i] << endl;
				//imshow(to_string(i), imagini[i]);
				similare[1] = imagini[i];


			}
			if (simili[i] == similiCop[index - 3])
			{
				//cout << endl;
				//cout << "img" << i << "___" << simili[i] << endl;
				//imshow(to_string(i), imagini[i]);
				similare[2] = imagini[i];

			}
		}
		
		//imshow("Patched Image", src);
		//waitKey();
	}

}*/

void select3Poze(Mat pozeSim[], float similaritati[],char fname[])
{
	Mat src = imread(fname, IMREAD_COLOR);
	Mat imagini[1000];
	Mat bucatiOriginal[1000];
	int indexOriginal = 0;
	float Voriginal[1000];
	int indexVOriginal;
	int index = 0;
	int N;
	cout << "M = ";
	cin >> N;
	float V[1000];
	int indexV;
	float simili[1000] = { 0.0 };
	gridImagine(src, bucatiOriginal, indexOriginal, N, Voriginal, indexVOriginal);
	testOpenImagesFldModif2(imagini, index);
	for (int i = 0; i < index; i++)
	{
		int index2 = 0;
		Mat bucati[1000] = {};
		gridImagine(imagini[i], bucati, index2, N, V, indexV);
		for (int j = 0; j < N * N; j++)
		{
			simili[i] = similitudini(Voriginal, V, N);
		}
	}
	float similiCop[100];
	for (int i = 0; i < index; i++)
	{
		similiCop[i] = simili[i];
	}
	float temp = 0.0f;
	for (int i = 0; i < index; i++)
	{
		for (int j = 0; j < index - i - 1; j++)
		{
			if (similiCop[j] > similiCop[j + 1])
			{
				temp = similiCop[j + 1];
				similiCop[j + 1] = similiCop[j];
				similiCop[j] = temp;
			}
		}
	}
	/*for (int i = 0; i < index; i++)
	{
		cout << "img" << i << "___" << simili[i]<<endl;
	}
	for (int i = 0; i < index; i++)
	{
		cout << "imgcop" << i << "___" << similiCop[i]<<endl;
	}*/
	for (int i = 0; i < index; i++)
	{
		if (simili[i] == similiCop[index - 1])
		{
			pozeSim[1] = imagini[i];
			similaritati[1] = simili[i];
		}
		if (simili[i] == similiCop[index - 2])
		{
			pozeSim[2] = imagini[i];
			similaritati[2] = simili[i];
		}
		if (simili[i] == similiCop[index - 3])
		{
			pozeSim[3] = imagini[i];
			similaritati[3] = simili[i];
		}
	}
	pozeSim[0] = src;
	similaritati[0] = 1;
}


Mat filtruGaussian(Mat src) {
		Mat srcR = Mat(src.rows, src.cols, CV_8UC1);
		Mat srcG = Mat(src.rows, src.cols, CV_8UC1);
		Mat srcB = Mat(src.rows, src.cols, CV_8UC1);
		for (int i = 0; i < src.rows; i++) 
		{
			for (int j = 0; j < src.cols; j++) 
			{
				Vec3b v3 = src.at<Vec3b>(i, j);
				uchar b = v3[0];
				uchar g = v3[1];
				uchar r = v3[2];
				srcR.at<uchar>(i, j) = r;
				srcG.at<uchar>(i, j) = g;
				srcB.at<uchar>(i, j) = b;
			}
		}

		Mat finalR = srcR.clone();
		Mat finalG = srcG.clone();
		Mat finalB = srcB.clone();
		Mat fin = Mat(src.rows, src.cols, CV_8UC3);
		int w = 5;
		double t = (double)getTickCount();

		float nucleuConvolutie[100][100];
		float row = (float)w / 6;
		int mid = w / 2;

		nucleuConvolutie[0][0] = nucleuConvolutie[4][4] = nucleuConvolutie[0][4] = nucleuConvolutie[4][0] = 0.0005;
		nucleuConvolutie[0][1] = nucleuConvolutie[0][3] = nucleuConvolutie[1][0] = nucleuConvolutie[1][4] = nucleuConvolutie[3][0] = nucleuConvolutie[3][4] = nucleuConvolutie[4][1] = nucleuConvolutie[4][3] = 0.0050;
		nucleuConvolutie[0][2] = nucleuConvolutie[2][0] = nucleuConvolutie[2][4] = nucleuConvolutie[4][2] = 0.0109;
		nucleuConvolutie[1][1] = nucleuConvolutie[1][3] = nucleuConvolutie[3][1] = nucleuConvolutie[3][3] = 0.0521;
		nucleuConvolutie[1][2] = nucleuConvolutie[2][1] = nucleuConvolutie[2][3] = nucleuConvolutie[3][2] = 0.1139;
		nucleuConvolutie[2][2] = 0.2487;

		for (int i = 0; i < w; i++) {
			for (int j = 0; j < w; j++) {
				//cout << nucleuConvolutie[i][j] << "    ";
			}
			//cout << endl;
		}

		for (int i = w / 2; i < src.rows - w / 2; i++) {
			for (int j = w / 2; j < src.cols - w / 2; j++) {

				float valoareR = 0;
				float valoareG = 0;
				float valoareB = 0;


				for (int u = 0; u < w; u++) {
					for (int v = 0; v < w; v++) {
						valoareR = valoareR + nucleuConvolutie[u][v] * srcR.at<uchar>(i + u - w / 2, j + v - w / 2);
						valoareG = valoareG + nucleuConvolutie[u][v] * srcG.at<uchar>(i + u - w / 2, j + v - w / 2);
						valoareB = valoareB + nucleuConvolutie[u][v] * srcB.at<uchar>(i + u - w / 2, j + v - w / 2);

					}
				}
				finalR.at<uchar>(i, j) = (char)(valoareR);
				finalG.at<uchar>(i, j) = (char)(valoareG);
				finalB.at<uchar>(i, j) = (char)(valoareB);
				fin.at<Vec3b>(i, j)[0] = finalB.at<uchar>(i, j);
				fin.at<Vec3b>(i, j)[1] = finalG.at<uchar>(i, j);
				fin.at<Vec3b>(i, j)[2] = finalR.at<uchar>(i, j);


			}
		}
		for (int i = w / 2; i < src.rows - w / 2; i++) {
			for (int j = w / 2; j < src.cols - w / 2; j++) {
				fin.at<Vec3b>(i, j)[0] = finalB.at<uchar>(i, j);
				fin.at<Vec3b>(i, j)[1] = finalG.at<uchar>(i, j);
				fin.at<Vec3b>(i, j)[2] = finalR.at<uchar>(i, j);
			}
		}
		return fin;
		//imshow("originalImage", src);
		//imshow("gaussianImage", fin);

		//waitKey();


}
void convHSI_histograma(Mat src, int hist[]) {

	Mat hsi = Mat(src.rows, src.cols, CV_8UC3);
	int r, g, b;
	double h, s, in;
	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			r = src.at<Vec3b>(i, j)[2];
			g = src.at<Vec3b>(i, j)[1];
			b = src.at<Vec3b>(i, j)[0];

			in = (b + g + r) * 1.00 / 3;

			int minimRGB = 0;
			minimRGB = min(r, min(b, g));
			if (in == 0)
				s = 0;
			else
				s = 1 - (minimRGB / in);

			if (s != 0)
			{
				h = 0.5 * ((r - g) + (r - b)) / sqrt(((r - g) * (r - g)) + ((r - b) * (g - b)));
				h = acos(h);

				if (b > g)
				{
					h = ((360 * PI) / 180.0) - h;
				}
				else
					h = h;
			}
			else {
				h = 0;
			}
			hsi.at<Vec3b>(i, j)[0] = ((h * 180) / PI) * 255 / 360;
			hsi.at<Vec3b>(i, j)[1] = s * 100;
			hsi.at<Vec3b>(i, j)[2] = in;
		}
	}
	hist[256] = { 0 };
	for (int i = 0; i < src.rows; i++)
	{
		for (int j = 0; j < src.cols; j++)
		{
			int pixel = hsi.at<Vec3b>(i, j)[0];
			hist[pixel]++;
		}
	}
	//return hist;
	//showHistogram("histograma", hist, 256, 200);
}

float correlationCoefficient(float X[], float Y[], int n)
{

	float sum_X = 0., sum_Y = 0., sum_XY = 0.;
	float squareSum_X = 0., squareSum_Y = 0.;

	for (int i = 0; i < n; i++)
	{
		sum_X += X[i];
		sum_Y += Y[i];
		sum_XY += X[i] * Y[i];
		squareSum_X += X[i] * X[i];
		squareSum_Y += Y[i] * Y[i];
	}
	/*
	cout << "sumX = " << sum_X << endl;
	cout << "sumY = " << sum_Y << endl;
	cout << "sumXY = " << sum_XY << endl;
	cout << "sumX2 = " << squareSum_X << endl;
	cout << "sumY2 = " << squareSum_Y << endl;
	*/
	float corr = (float)(n * sum_XY - (sum_X * sum_Y)) / sqrt((n * squareSum_X - (sum_X * sum_X)) * (n * squareSum_Y - (sum_Y * sum_Y)));
	return corr;
}

void proiect2()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat pozeSim[4];
		float similaritati[4];
		select3Poze(pozeSim, similaritati, fname);
		imshow("originala", pozeSim[0]);
		/*imshow("similara1", pozeSim[1]);
		imshow("similara2", pozeSim[2]);
		imshow("similara3", pozeSim[3]);*/
		cout << "similitudine DB ==> " << similaritati[1]<<", "<< similaritati[2]<<", "<< similaritati[3] << endl;
		Mat blurOrig = Mat(pozeSim[0].rows, pozeSim[0].cols, CV_8UC3);
		Mat blur1 = Mat(pozeSim[1].rows, pozeSim[1].cols, CV_8UC3);
		Mat blur2 = Mat(pozeSim[2].rows, pozeSim[2].cols, CV_8UC3);
		Mat blur3 = Mat(pozeSim[3].rows, pozeSim[3].cols, CV_8UC3);
		blurOrig = filtruGaussian(pozeSim[0]);
		blur1 = filtruGaussian(pozeSim[1]);
		blur2 = filtruGaussian(pozeSim[2]);
		blur3 = filtruGaussian(pozeSim[3]);
		//imshow("originala", blurOrig);
		//imshow("similara1", blur1);
		//imshow("similara2", blur2);
		//imshow("similara3", blur3);
		int histOrig[256] = { 0 };
		int hist1[256] = { 0 };
		int hist2[256] = { 0 };
		int hist3[256] = { 0 };

		convHSI_histograma(blurOrig, histOrig);
		convHSI_histograma(blur1, hist1);
		convHSI_histograma(blur2, hist2);
		convHSI_histograma(blur3, hist3);
		showHistogram("Histograma", histOrig, 256, 200);
		showHistogram("Histograma1", hist1, 256, 200);
		showHistogram("Histograma2", hist2, 256, 200);
		showHistogram("Histogram3", hist3, 256, 200);
		
		int n = 256;
		float histOrig2[256] = { 0. };
		float hist12[256] = { 0. };
		float hist22[256] = { 0. };
		float hist32[256] = { 0. };
		for (int i = 0; i < n; i++) {
			//cout << (float)histOrig[i]/(blurOrig.rows*blurOrig.cols)<<" ";
			histOrig2[i] = (float)histOrig[i] / (blurOrig.rows * blurOrig.cols);
			hist12[i] = (float)hist1[i] / (blur1.rows * blur1.cols);
			hist22[i] = (float)hist2[i] / (blur2.rows * blur2.cols);
			hist32[i] = (float)hist3[i] / (blur3.rows * blur3.cols);

		}
		//cout << "coef1" << endl;
		float coef1 = correlationCoefficient(histOrig2, hist12, n);
		//cout << "coef2" << endl;
		float coef2 = correlationCoefficient(histOrig2, hist22, n);
		//cout << "coef3" << endl;
		float coef3 = correlationCoefficient(histOrig2, hist32, n);
		
		float coeficienti[4] = {0.0};
		float coefOrd[4] = {0.0};
		float temp = 0.0;
		coeficienti[0] = coef1;
		coeficienti[1] = coef2;
		coeficienti[2] = coef3;
		coefOrd[0] = coef1;
		coefOrd[1] = coef2;
		coefOrd[2] = coef3;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3 - i - 1; j++)
			{
				if (coefOrd[j] < coefOrd[j + 1])
				{
					temp = coefOrd[j + 1];
					coefOrd[j + 1] = coefOrd[j];
					coefOrd[j] = temp;
				}
			}
		}
		
		cout << "Corelare Pearson ==> " << coefOrd[0] << ", " << coefOrd[1] << ", " << coefOrd[2] << endl;

		for (int i = 0; i < 4; i++)
		{
			if (coeficienti[i] == coefOrd[0])
			{
				imshow("similara1", pozeSim[i+1]);
			}
			if (coeficienti[i] == coefOrd[1])
			{
				imshow("similara2", pozeSim[i+1]);
			}
			if (coeficienti[i] == coefOrd[2])
			{
				imshow("similara3", pozeSim[i+1]);
			}
		}
		waitKey();
	}

}

int main()
{
	proiect2();
	return 0;
}