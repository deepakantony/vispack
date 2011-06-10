#include <imagefile.h>
#include "edgedetector.h"
#include "util.h"
#include "houghtransform.h"


int main(int argc, char* argv[]) {
	if(argc < 11) {
		cout << "Usage: <program> <image> <median-window-size> <median-iterations> <sigma> <threshold> <min-radius> <max-radius> <radius-step> <grey-value> <ht-threshold> <output_filename>" << endl;
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
	
	writeImage(op, "edge.tif");
	
	float minradius;
	sscanf(argv[6], "%f", &minradius);
	
	float maxradius;
	sscanf(argv[7], "%f", &maxradius);
	
	float radiusstep;
	sscanf(argv[8], "%f", &radiusstep);
	
	float greyValue;
	sscanf(argv[9], "%f", &greyValue);
	updateGreyVal (greyValue);
	
	float ht_threshold;
	sscanf(argv[10], "%f", &ht_threshold);
	
	VISImage<float> opHT = circleHT(op, med, minradius, maxradius, radiusstep, ht_threshold);
	
	if(argc > 11)
		writeImage(opHT, argv[11]);
	else writeImage(opHT, "htcircleswithradius.tif");
}
