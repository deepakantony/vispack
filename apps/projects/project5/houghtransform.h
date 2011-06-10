#ifndef _HOUGHTRANSFORM_H
#define _HOUGHTRANSFORM_H

#include <image.h>
#include <array.h>
#include <matrix.h>
#include <assert.h>
#include "util.h"

#define DEFAULT_RADIUS -1.f
#define DEFAULT_THRESHOLD 0.7f
#define DEFAULT_TOLERANCE 2.0f

// update the gray value
void updateGreyVal(float newGV);

// sets up accumulator array
// x, y & z are the size of the accumulator in each dimension
VISArray< VISImage<int> > get3DAccumulator(int x, int y, int z);

// set up the max hit for the accumulator
VISArray< VISImage<int> > getMaxHitForAcc(VISArray< VISImage<int> > &acc, VISVector accSize);

// updates accumulator array
// acc is the accumulator, u & v are the edge point in the image, accSize gives the size in each direction
// rStep gives the step value for the radius
void updateCircleAcc(VISArray< VISImage<int> > &acc, int u, int v, VISVector accSize, float rStep, float tolerance);

// find the circles
VISArray<VISVector> findCircles(VISArray< VISImage<int> > &acc, VISVector accSize, float rStep, float threshold, int &nCircles);

// draw the circle in the image
void drawCircle(VISImage<float> &I, int x0, int y0, float radius);

// Get the final image with the circles 
VISImage<float> finalImage(VISImage<float> I, VISArray<VISVector> circles, int nCircles);

// Maximum radius is half the diagonal of the image
float maxRadius(int w, int h);

// Returns the circles found through hough transform
VISArray<VISVector> houghTransformCircles(VISImage<float> I, float radius, float threshold, int &nCircles);

// hough transform for circles with radius
VISImage<float> circleHTwithRad(VISImage<float> I, VISImage<float> origImage, float radius, float threshold = DEFAULT_THRESHOLD);

// hough transform for circles without radius
VISImage<float> circleHT(VISImage<float> I, VISImage<float> origImage, float minRadius, float maxRadius, float radiusStep, float threshold = DEFAULT_THRESHOLD);

#endif /* HOUGHTRANSFORM_H */

