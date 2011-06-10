#include <imagefile.h>
#include "edgedetector.h"
#include "util.h"


int main(int argc, char* argv[]) {
	if(argc < 6) {
		cout << "Usage: <program> <image> <median-window-size> <median-iterations> <sigma> <threshold> <output_filename>" << endl;
		return -1;
	}
	
	VISIm im;
	VISImageFile imgFile;
	VISImage<float> I;
	if((im = imgFile.read(argv[1])).isValid()) {
		I = VISImage<float>(im);
	}
	
	float sigma;
	sscanf(argv[4], "%f", &sigma);
	
	float threshold = 0.01f;
	sscanf(argv[5], "%f", &threshold);
	
	// Pass it through a median filter
	float its;
	float winsize;
	sscanf(argv[2], "%f", &winsize);
	sscanf(argv[3], "%f", &its);
	VISImage<float> med = I;
	for(int i = 0; i < its; i++)
		med = med.median((int)winsize);
	writeImage(med, "median.tif");
	
	VISImage<float> op = marrHildreth (med, sigma, threshold);
	
	if(argc > 6)
		writeImage(op, argv[6]);
	else writeImage(op, "edge.tif");
}
