#ifndef _UTIL_H
#define _UTIL_H

#include <array.h>
#include <param.h>
#include <imagefile.h>
#include <matrix.h>
#include <image.h>
#include <cmath>
#include <string>

// Checks if the array is empty
bool isEmpty(VISArray<bool> A);

// checks if the string ends with the string in end
bool endswith(string str, string end);

// returns the euclidean distance between p1 and p2
float euclideanDistance(VISVector p1, VISVector p2);

// uses euclidean distance to calculate the weight function
float weightFunction(VISVector point, VISVector center);

// Update minimum (x,y) from the given vector "c"
void updateMinsMaxs(VISVector c, int &minx, int &miny, int &maxx, int &maxy);

// A helper function that multiplies i'th element in vec1 with i'th ele in vec2.
// returns the resultant vector
VISVector mul(VISVector vec1, VISVector vec2);

// prints out the array
template<class T>
void displayArray(VISArray<T> array);

// A helper function which reads, from a ParameterFile object, a
// matrix of the form (name rows cols ....data in row major order...)
// Returns an invalid matrix if it fails (prints something)
VISMatrix readMatrix(VPF::ParameterFile pf, const char* name);

// Given an image in float array this function returns a VISImage equivalent object
VISImage<float> getVISImage(float *I, int w, int h);

// Given a VISImage object this function returns a float array equivalent
float* getFloatArray(VISImage<float> I);

void writeImage(VISImage<float> I, char* filename);

#endif /* UTIL_H */

