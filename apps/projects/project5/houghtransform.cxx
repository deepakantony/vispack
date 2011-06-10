#include "houghtransform.h"

float greyValue = 255.f;

// sets up accumulator array
// x, y & z are the size of the accumulator in each dimension
VISArray< VISImage<int> > get3DAccumulator(int x, int y, int z) {
	VISArray< VISImage<int> > acc;
	for(int i = 0; i < z; i++) {
		VISImage<int> img(x,y);
		img = 0;
		acc.poke(i) = img;
	}
	
	return acc;
}

// set up the max hit for the accumulator
VISArray< VISImage<int> > getMaxHitForAcc(VISArray< VISImage<int> > &acc, VISVector accSize) {
}

// updates accumulator array
// acc is the accumulator, u & v are the edge point in the image, accSize gives the size in each direction
// rStep gives the step value for the radius
void updateCircleAcc(VISArray< VISImage<int> > &acc, int u, int v, VISVector accSize, float rStep, float tolerance) {
/*	for(int x = 0; x < accSize.peek(0); x++)
		for(int y = 0; y < accSize.peek(1); y++) {
		float lhs = (u-x)*(u-x) + (v-y)*(v-y);
		if( lhs < (r+tolerance) && lhs > (r-tolerance) )
			acc.poke(z).poke(x,y) += 1;
		}*/
	
	float r = rStep;
	for(int z = 0; z < accSize.peek(2); z++, r += rStep) {		
		// Use bersenham algorithm to increment values in accumulator
		int f = 1 - r;
		int ddF_x = 1;
		int ddF_y = -2 * r;
		int x = 0;
		int y = r;

		if((v+r) < accSize.peek(1))
			acc.poke(z).poke(u, v + r) += 1;
		if((v-r) >= 0)
			acc.poke(z).poke(u, v - r) += 1;
		if((u+r) < accSize.peek(0))
			acc.poke(z).poke(u + r, v) += 1;
		if((u-r) >= 0) 
			acc.poke(z).poke(u - r, v) += 1;

		while(x < y) {
			assert(ddF_x == 2 * x + 1);
			assert(ddF_y == -2 * y);
			assert(f == x*x + y*y - r*r + 2*x - y + 1);
			if(f >= 0) {
				y--;
				ddF_y += 2;
				f += ddF_y;
			}
			x++;
			ddF_x += 2;
			f += ddF_x;
			
			if((u+x) < accSize.poke(0)) {
				if((v+y) < accSize.poke(1))
					acc.poke(z).poke(u + x, v + y) += 1;
				if((v-y) >= 0)
					acc.poke(z).poke(u + x, v - y) += 1;
			}
			
			if((u+y) < accSize.poke(0)) {
				if((v+x) < accSize.poke(1))
					acc.poke(z).poke(u + y, v + x) += 1;
				if((v-x) >= 0)
					acc.poke(z).poke(u + y, v - x) += 1;
			}
			
			if((u-x) >= 0) {
				if((v+y) < accSize.poke(1))
					acc.poke(z).poke(u - x, v + y) += 1;
				if((v-y) >= 0)
					acc.poke(z).poke(u - x, v - y) += 1;
			}
			
			if((u-y) >= 0) {
				if((v+x) < accSize.poke(1))
					acc.poke(z).poke(u - y, v + x) += 1;
				if((v-x) >= 0)
					acc.poke(z).poke(u - y, v - x) += 1;
			}
		}
	}
}

// find the circles
VISArray<VISVector> findCircles(VISArray< VISImage<int> > &acc, VISVector accSize, float rStep, float threshold, int &nCircles) {
	int max = 0; // stores the max value that will be used to threshold the accumulator
	
	// First I apply the method of suppressing non-maximum local values
	// i.e. check every cell if it has a higher value than its neighbours
	// if true the value is left alone otherwise it is changed to 0
	for(int x = 0; x < (int)accSize.peek(0); x++)
		for(int y = 0; y < (int)accSize.peek(1); y++) {
			for(int z = 0; z < (int)accSize.peek(2); z++) {
				bool isMax = true;
				int curVal = acc.peek(z).peek(x,y);
				
				for(int xx = (x == 0?0:(x-1)); (xx < (x+1)) && (xx < accSize.peek(0)); xx++)
					for(int yy = (y == 0?0:(y-1)); (yy < (y+1)) && (yy < accSize.peek(1)); yy++)
						for(int zz = (z == 0?0:(z-1)); (zz < (z+1)) && (zz < accSize.peek(2)); zz++)
							if(acc.peek(zz).peek(xx,yy) >= curVal && !(xx == x && yy == y && zz == z))
								isMax = false;
				
				if(!isMax)
					acc.poke(z).poke(x,y) = 0;
				else if (curVal > max)
					max = curVal;
			}
		}

	if(threshold > 1 || threshold < 0)
		threshold = DEFAULT_THRESHOLD;
	int T = max*threshold;

	VISArray<VISVector> ret;
	int i = 0;
	nCircles = 0;
	for(int x = 0; x < accSize.peek(0); x++)
		for(int y = 0; y < accSize.peek(1); y++) {
			float r = rStep;
			for(int z = 0; z < (int)accSize.peek(2); z++, r += rStep)
				if(acc.peek(z).peek(x,y) > T) {
					VISVector c(x,y,r);
					ret.poke(i) = c;
					i++;
					nCircles++;
				}
		}

	return ret;
}

// draw the circle in the image
// bresenham circle algorithm taken from wiki
void drawCircle(VISImage<float> &I, int x0, int y0, float radius) {
	int w = I.width();
	int h = I.height();

	int f = 1 - radius;
	int ddF_x = 1;
	int ddF_y = -2 * radius;
	int x = 0;
	int y = radius;

	if( (y0+radius) < h )
		I.poke(x0, y0 + radius) = greyValue;
	if( (y0-radius) >= 0 )
		I.poke(x0, y0 - radius) = greyValue;
	if( (x0+radius) < w )
		I.poke(x0 + radius, y0) = greyValue;
	if( (x0-radius) >= 0 )
		I.poke(x0 - radius, y0) = greyValue;

	while(x < y)
	{
		assert(ddF_x == 2 * x + 1);
		assert(ddF_y == -2 * y);
		assert(f == x*x + y*y - radius*radius + 2*x - y + 1);
		if(f >= 0) 
		{
			y--;
			ddF_y += 2;
			f += ddF_y;
		}
		x++;
		ddF_x += 2;
		f += ddF_x;
		
		if((x0+x) < w) {
			if((y0+y) < h)
				I.poke(x0 + x, y0 + y) = greyValue;
			if((y0-y) >= 0)
				I.poke(x0 + x, y0 - y) = greyValue;
		}
		
		if((x0+y) < w) {
			if((y0+x) < h)
				I.poke(x0 + y, y0 + x) = greyValue;
			if((y0-x) >= 0)
				I.poke(x0 + y, y0 - x) = greyValue;
		}
		
		if((x0-x) >= 0) {
			if((y0+y) < h)
				I.poke(x0 - x, y0 + y) = greyValue;
			if((y0-y) >= 0)
				I.poke(x0 - x, y0 - y) = greyValue;
		}
		
		if((x0-y) >= 0) {
			if((y0+x) < h)
				I.poke(x0 - y, y0 + x) = greyValue;
			if((y0-x) >= 0)
				I.poke(x0 - y, y0 - x) = greyValue;
		}
	}
}

// Get the final image with the circles 
VISImage<float> finalImage(VISImage<float> I, VISArray<VISVector> circles, int nCircles){
	VISImage<float> ret(I);
	
	for(int i = 0; i < nCircles; i++) {
		drawCircle (ret, circles.peek(i).peek(0), circles.peek(i).peek(1), circles.peek(i).peek(2));
	}
	
	return ret;
}

// Maximum radius is half the diagonal of the image
float maxRadius(int w, int h) {
	return 0.5f*sqrt(w*w+h*h);
}

// Returns the circles found through hough transform
VISArray<VISVector> houghTransformCircles(VISImage<float> I, float radius, float threshold, int &nCircles) {
	int w = I.width();
	int h = I.height();

	VISVector accSize(3);
	accSize.poke(0) = w;
	accSize.poke(1) = h;
	accSize.poke(2) = 1;
	cout << "Accumulator size initialized!" << endl;
	
	float rStep = radius;
	float tolerance = DEFAULT_TOLERANCE;
	
	// set up the accumulator
	VISArray< VISImage<int> > acc = get3DAccumulator (w,h,1);
	cout << "Accumulator set up!" << endl;
	
	// For each pixel in the image update the accumulator
	for(int u = 0; u < w; u++)
		for(int v = 0; v < h; v++)
			if(I.peek(u,v) > 0)
				updateCircleAcc (acc, u, v, accSize, rStep, tolerance);
	cout << "Accumulator data is updated." << endl;
	
	// blur the accumulator space
	//acc.poke(0) = acc.peek(0).convolve(gauss_kernel(1.f, 2));
	writeImage(acc.peek(0), "htacc.tif");
	
	// Find the circles in the accumulator
	nCircles = 0;
	VISArray<VISVector> circles = findCircles (acc, accSize, rStep, threshold, nCircles);
	cout << "Number of circles found: " << nCircles << endl;
	cout << "Found the circles from accumulator." << endl;
	
	return circles;
}

// update the gray value
void updateGreyVal(float newGV) {
	greyValue = newGV;
}

// hough transform for circles with radius
VISImage<float> circleHTwithRad(VISImage<float> I, VISImage<float> origImage, float radius, float threshold) {
	int nCircles = 0;	
	VISArray<VISVector> circles = houghTransformCircles(I, radius, threshold, nCircles);
	
	// Use the circles found to update the image
	return finalImage (origImage, circles, nCircles);
}

// hough transform for circles without radius
VISImage<float> circleHT(VISImage<float> I, VISImage<float> origImage, float minRadius, float maxRadius, float radiusStep, float threshold) {
	cout << "Min radius: " << minRadius << endl;
	cout << "Max radius: " << maxRadius << endl;
	cout << "Radius step: " << radiusStep << endl;
	VISImage<float> oi = origImage;
	if(radiusStep <= 0) radiusStep = 1.f;
	for(float i = minRadius; i <= maxRadius; i+=radiusStep) {
		cout << "Performing HT for radius: " << i << endl;
		int nCircles = 0;
		VISArray<VISVector> circles = houghTransformCircles(I, i, threshold, nCircles);
		
		// Use the circles found to update the image
		oi = finalImage (oi, circles, nCircles);
	}
	
	return oi;
}
