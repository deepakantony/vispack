#include <fstream>
#include <iostream>
#include <string>
#include <time.h>

#include <image.h>
#include <imageRGBA.h>
#include <imagefile.h>
#include <rgba.h>

using namespace std;

template <class T>
float *histogram(VISImage<T> I, int n, T min, T max) {
	float *H = (float *)malloc(sizeof(float)*n);
	
	for(int i = 0; i < n; i++)
		H[i] = 0.0f;
	
	T mulfact = (n-1)/(max-min);
	int h = I.height();
	int w = I.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			T intensity = I.peek(i,j);
			if( intensity >= min && intensity <= max)
				H[(int)((intensity-min)*mulfact)] += 1.0f;
		}
	
	return H;
}

template <class T>
VISImage<T> histeq(VISImage<T> I, int n, T min, T max, float alpha = 1.0f) {
	if(alpha < 0.0f) alpha = 0.0f;
	else if(alpha > 1.0f) alpha = 1.0f;
	
	// get the histogram
	float *H = histogram (I, n, min, max);

	// get the cumulative histogram values
	for(int i = 1; i < n; i++)
		H[i] = H[i]+H[i-1];
	
	int h = I.height();
	int w = I.width();
	int size = h*w;
	VISImage<T> modI = I.createToSize();
	modI = (T)0;
	
	T mulfact = (n-1)/(max-min);
	
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			T intensity = I.peek(i,j);
			T newIntensity = alpha * H[(int)((intensity-min)*mulfact)]*(n-1)/size + (1.0f-alpha)*intensity;
			
			if(newIntensity > max)
				newIntensity = max;
			if(newIntensity < min)
				newIntensity = min;
			
			modI.poke(i,j) = newIntensity;
		}
	
	return modI;
}

template <class T>
VISImage<T> AHE(VISImage<T> I, int n, T min, T max, int windowsize) {	
	int h = I.height();
	int w = I.width();
	int size = h*w;
	VISImage<T> modI = I.createToSize();
	modI = (T)0;
	
	T mulfact = (n-1)/(max-min);
	int wsSQR = windowsize*windowsize;
	T mulfact2 = ((T)(n-1))/((T)wsSQR);
	
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			T intensity = I.peek(i,j);
			
			// Initialize histogram to 0
			int *H = (int*)malloc(n*sizeof(int));
			for(int ii = 0; ii < n; ii++)
				H[ii] = 0;
			
			for(int c = i-windowsize/2; c < (i+windowsize/2); c++)
				for(int r = j-windowsize/2; r < (j+windowsize/2); r++) {
					if(c >= 0 && r >= 0 && c < w && r < h) {
						T l = I.peek(c, r);
						H[(int)((l-min)*mulfact)] += 1;
					}
					else {
						H[(int)((0-min)*mulfact)] += 1;
					}
				}
			
			int sum = 0;
			for(int ii = 0; ii <= ((int)intensity); ii++)
				sum += H[ii];
			
			T l = sum*mulfact2;
			if(l > max) l = max;
			else if(l < min) l = min;
			modI.poke(i,j) = l;
			
			free(H);
		}
	
	return modI;
}

VISImage<float> AHE(VISImage<float> I, int n, int windowsize = -1) {
	return AHE(I, n, 0.0f, 255.0f, windowsize);
}

template <class T>
VISImage<T> CLAHE(VISImage<T> I, int n, T min, T max, int cliplevel, int windowsize) {	
	int h = I.height();
	int w = I.width();
	int size = h*w;
	VISImage<T> modI = I.createToSize();
	modI = (T)0;
	
	T mulfact = (n-1)/(max-min);
	int wsSQR = windowsize*windowsize;
	T mulfact2 = ((T)(n-1))/((T)wsSQR);
	
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			T intensity = I.peek(i,j);
			
			// Initialize histogram to 0
			int *H = (int*)malloc(n*sizeof(int));
			for(int ii = 0; ii < n; ii++)
				H[ii] = 0;
			
			for(int c = i-windowsize/2; c < (i+windowsize/2); c++)
				for(int r = j-windowsize/2; r < (j+windowsize/2); r++) {
					if(c >= 0 && r >= 0 && c < w && r < h) {
						T l = I.peek(c, r);
						H[(int)((l-min)*mulfact)] += 1;
					}
					else {
						H[(int)((0-min)*mulfact)] += 1;
					}
				}
			
			int clippedpixels = 0;
			for(int ii = 0; ii < n; ii++) {
				if(H[ii] > cliplevel) {
					clippedpixels += (H[ii]-cliplevel);
					H[ii] = cliplevel;
				}
			}
			
			int idx = 0;
			while(clippedpixels > 0) {
				H[idx]++;
				clippedpixels--;
				
				idx = (idx+1)%n;
			}
			
			int sum = 0;
			for(int ii = 0; ii < ((int)intensity); ii++)
				sum += H[ii];
			T l = sum*mulfact2;
			if(l > max) l = max;
			else if(l < min) l = min;
			
			modI.poke(i,j) = l;
			
			free(H);
		}
	
	return modI;
}

VISImage<float> CLAHE(VISImage<float> I, int n, int cliplevel, int windowsize = -1) {
	return CLAHE(I, n, 0.0f, 255.0f, cliplevel, windowsize);
}

void menuHistEQ(VISImage<float> I, char *filename) {
	cout << "Global Histogram Equalization" << endl;
	cout << "Please enter the parameters for histogram equalization" << endl;
	
	cout << "Number of bins: ";
	int n;
	cin >> n;
	cout << "Minimum intensity value: ";
	float min;
	cin >> min;
	cout << "Maximum intensity value: ";
	float max;
	cin >> max;
	cout << "Alpha (blending value between 0 and 1), anything outside the limit will be treated as 1: ";
	float alpha;
	cin >> alpha;
	
	if(alpha < 0.0f || alpha > 1.0f)
		alpha = 1.0f;
	
	
	float *origHist = histogram (I, n, min, max);
	VISImage<float> newIf(histeq(I, n, min, max, alpha));
	float *newHist = histogram (newIf, n, min, max);
	
	// convert float to byte
	VISImage<byte> newI(newIf);
	
	// fn = filename
	string imgfn(filename);
	imgfn.append("_histeq.tif");
	string histfn(filename);
	histfn.append("_hist.txt");
	string histeqfn(filename);
	histeqfn.append("_histeq.txt");
	
	// write the image file
	VISImageFile imFile;
	imFile.write_tiff(newI, imgfn.c_str());
	
	// write the histograms into files
	ofstream output;
	output.open(histfn.c_str());
	for(int i = 0; i < (n-1); i++)
		output << origHist[i] << "\n";
	output << origHist[n-1] << endl;
	output.close();
	
	output.open(histeqfn.c_str());
	for(int i = 0; i < (n-1); i++)
		output << newHist[i] << "\n";
	output << newHist[n-1] << endl;
	output.close();
	
	cout << "Equalized image written into " << imgfn << endl;
	cout << "Histograms written into " << histfn << " (original image) and " 
		<< histeqfn << " (modified image)." << endl;

}

template <class T>
void menuAHE(VISImage<T> I, char *filename) {
	cout << "Adaptive Histogram Equalization" << endl;
	cout << "Please enter the parameters for histogram equalization" << endl;
	
	cout << "Number of bins: ";
	int n;
	cin >> n;
	cout << "Window size: ";
	float windSize;
	cin >> windSize;
	
	cout << "Please wait.. This might take a while! ...." << endl;
	
	if(windSize < 0)
		windSize = I.width()/8;
	
	float *origHist = histogram (I, n, 0.0f, 255.0f);
	time_t start,end;
	time (&start);
	VISImage<float> newIf = AHE(I, n, windSize);
	time (&end);
	double dif = difftime (end,start);
	cout << "It took " << dif << " seconds to execute AHE" << endl;
	float *newHist = histogram (newIf, n, 0.0f, 255.0f);
	
	// convert float to byte
	VISImage<byte> newI(newIf);
	
	// fn = filename
	string imgfn(filename);
	imgfn.append("_AHE.tif");
	string histfn(filename);
	histfn.append("_hist.txt");
	string histeqfn(filename);
	histeqfn.append("_AHE.txt");
	
	// write the image file
	VISImageFile imFile;
	imFile.write_tiff(newI, imgfn.c_str());
	
	// write the histograms into files
	ofstream output;
	output.open(histfn.c_str());
	for(int i = 0; i < (n-1); i++)
		output << origHist[i] << "\n";
	output << origHist[n-1] << endl;
	output.close();
	
	output.open(histeqfn.c_str());
	for(int i = 0; i < (n-1); i++)
		output << newHist[i] << "\n";
	output << newHist[n-1] << endl;
	output.close();
	
	cout << "Equalized image written into " << imgfn << endl;
	cout << "Histograms written into " << histfn << " (original image) and " 
		<< histeqfn << " (modified image)." << endl;
}

template <class T>
void menuCLAHE(VISImage<T> I, char *filename) {
	cout << "Clipped Local Adaptive Histogram Equalization" << endl;
	cout << "Please enter the parameters for histogram equalization" << endl;
	
	cout << "Number of bins: ";
	int n;
	cin >> n;
	
	cout << "Enter the clipping limit between 0 and 1: ";
	float cliplim;
	cin >> cliplim;
	
	cout << "Window size: ";
	float windSize;
	cin >> windSize;
	
	cout << "Please wait.. This might take a while! ...." << endl;
	
	if(windSize < 0)
		windSize = I.width()/8;
	
	float *origHist = histogram (I, n, 0.0f, 255.0f);
	time_t start,end;
	time (&start);
	VISImage<float> newIf = CLAHE(I, n, cliplim*windSize*windSize, windSize);
	time (&end);
	double dif = difftime (end,start);
	cout << "It took " << dif << " seconds to execute CLAHE" << endl;
	float *newHist = histogram (newIf, n, 0.0f, 255.0f);
	
	// convert float to byte
	VISImage<byte> newI(newIf);
	
	// fn = filename
	string imgfn(filename);
	imgfn.append("_CLAHE.tif");
	string histfn(filename);
	histfn.append("_hist.txt");
	string histeqfn(filename);
	histeqfn.append("_CLAHE.txt");
	
	// write the image file
	VISImageFile imFile;
	imFile.write_tiff(newI, imgfn.c_str());
	
	// write the histograms into files
	ofstream output;
	output.open(histfn.c_str());
	for(int i = 0; i < (n-1); i++)
		output << origHist[i] << "\n";
	output << origHist[n-1] << endl;
	output.close();
	
	output.open(histeqfn.c_str());
	for(int i = 0; i < (n-1); i++)
		output << newHist[i] << "\n";
	output << newHist[n-1] << endl;
	output.close();
	
	cout << "Equalized image written into " << imgfn << endl;
	cout << "Histograms written into " << histfn << " (original image) and " 
		<< histeqfn << " (modified image)." << endl;
}

int main(int argc, char** argv) {
	VISIm im;
	VISImageFile imgFile;
	VISImage<float> I;
	if(argc > 1 && (im = imgFile.read(argv[1])).isValid()) {
		I = VISImage<float>(im);
	}
	else {
		cout << "Image file could not be read or no file was specified!" << endl;
		return -1;
	}
	cout << "Image's size: " << I.height() << "x" << I.width() << endl;
	cout << "Image's max intensity is " << I.max() << " and min intensity is " << I.min() << endl;
	
	while(true) {
		cout << "Here are your options... Please enter the corresponding number." << endl;
		
		cout << "1. Histogram Equalization" << endl;
		cout << "2. Adaptive histogram equalization" << endl;
		cout << "3. Clipped local adaptive histogram equalization" << endl;
		cout << "4. Exit" << endl;
		
		int option;
		cin >> option;
		
		switch(option) {
			case 1: menuHistEQ(I, argv[1]); break;
			case 2: menuAHE(I, argv[1]); break;
			case 3: menuCLAHE(I, argv[1]); break;
			case 4: return 0;
			default: 
				cout << "Wrong option entered!" << endl;
				break;
		}
	}
	return 0;
}
