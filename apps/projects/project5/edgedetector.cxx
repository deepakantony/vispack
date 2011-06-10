#include "edgedetector.h"

// Marr-Hildreth edge detector
VISImage<float> marrHildreth(VISImage<float> I, float sigma, float threshold) {	
	// smoothen with a gaussian filter
	VISImage<float> I_gauss = I.gauss(sigma);
	
	// find the laplacian of the image
	VISImage<float> lap(3,3);
	lap = -1.f;
	lap.poke(1,1) = 8.f;
	VISImage<float> I_lap = I_gauss.convolve(lap);
	
	cout << "laplacian of I, max: " << I_lap.max() << ", min: " << I_lap.min() << endl;
	
	// another way to find laplacian!
	//I_gauss.convolve(lap);
	//VISImage<float> I_dx = I_gauss.dx(2);
	//VISImage<float> I_dy = I_gauss.dy(2);
	//VISImage<float> I_lap = I_dx+I_dy;
	
	// Find the zero-crossings
	return zeroCrossings(I_lap, threshold*I_lap.max());
	//return I_lap
}

// Finds the zero crossings for the given image with a threshold on the gradient values
VISImage<float> zeroCrossings(VISImage<float> I, float threshold) {
	int w = I.width();
	int h = I.height();
	
	VISImage<float> retI(w,h);
	retI = 0.f;
	
	for(int i = 1; i < w-1; i++)
		for(int j = 1; j < h-1; j++) {
			if(I.peek(i,j) > threshold)
				if(I.peek(i-1, j) < 0 || I.peek(i+1, j) < 0 || I.peek(i, j-1) < 0 || I.peek(i, j+1) < 0)
					retI.poke(i,j) = 127.f;
		}
	
	return retI;
}

// Finds the zero crossings in the given image
// Marks it black on a white background
VISImage<float> zeroCrossings(VISImage<float> I) {
	int w = I.width();
	int h = I.height();
	
	VISImage<float> retI(w,h);
	retI = 0.f;
	
	for(int i = 1; i < w-1; i++)
		for(int j = 1; j < h-1; j++) {
			// First pair: left and right
			if( (I.peek(i-1, j)*I.peek(i+1, j)) < 0 )
				retI.poke(i,j) = 255.f;
			// Second pair: above and below
			else if( (I.peek(i, j-1)*I.peek(i, j+1)) < 0 )
				retI.poke(i,j) = 255.f;
			// Third pair: left above and right below
			//else if( (I.peek(i-1, j-1)*I.peek(i+1, j+1)) < 0 )
				//retI.poke(i,j) = 255.f;
			// Fourth pair: right above and left below
			//if( (I.peek(i+1, j-1)*I.peek(i-1, j+1)) < 0 )
				//retI.poke(i,j) = 255.f;
		}
	
	return retI;
}
