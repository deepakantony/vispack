#include "util.h"

// Checks if the array is empty
bool isEmpty(VISArray<bool> A) {
	for(int i = 0; i < A.n(); i++)
		if(A.peek(i))
			return false;
	return true;
}

// checks if the string ends with the string in end
bool endswith(string str, string end) {
	string::reverse_iterator ritend = end.rbegin();
	string::reverse_iterator ritstr = str.rbegin();
	while(ritend < end.rend() && ritstr < str.rend()) {
		if(*ritend != *ritstr)
			return false;
		
		ritend++;
		ritstr++;
	}
	return true;
}

// returns the euclidean distance between p1 and p2
float euclideanDistance(VISVector p1, VISVector p2) {
	float x = (p1.peek(0) - p2.peek(0))*(p1.peek(0) - p2.peek(0));
	float y = (p1.peek(1) - p2.peek(1))*(p1.peek(1) - p2.peek(1));
	return sqrt(x+y);
}

// uses euclidean distance to calculate the weight function
float weightFunction(VISVector point, VISVector center) {
	float distance = euclideanDistance(point, center);
	VISVector origin(2);
	origin = 0;
	float maxdistance = euclideanDistance(center, origin);
	
	return maxdistance-distance;
}

// Update minimum (x,y) from the given vector "c"
void updateMinsMaxs(VISVector c, int &minx, int &miny, int &maxx, int &maxy) {
	if(c.peek(0) < minx)
		minx = c.peek(0);
	if(c.peek(0) > maxx)
		maxx = c.peek(0);
	
	if(c.peek(1) < miny)
		miny = c.peek(1);
	if(c.peek(1) > maxy)
		maxy = c.peek(1);
}

// A helper function that multiplies i'th element in vec1 with i'th ele in vec2.
// returns the resultant vector
VISVector mul(VISVector vec1, VISVector vec2) {
	if(vec1.n() != vec2.n()) {
		cout << "Cannot multiply: the size of the vectors doesn't match" << endl;
		return VISVector();
	}
	
	VISVector result(vec1.n());
	for(int i = 0; i < vec1.n(); i++)
		result.poke(i) = vec1.peek(i)*vec2.peek(i);
	
	return result;
}

// prints out the array
template<class T>
void displayArray(VISArray<T> array) {
	for(int i = 0; i < array.n(); i++)
		cout << array[i] << " ";
	cout << endl;
}

// A helper function which reads, from a ParameterFile object, a
// matrix of the form (name rows cols ....data in row major order...)
// Returns an invalid matrix if it fails (prints something)
VISMatrix readMatrix(VPF::ParameterFile pf, const char* name)
{
  int r, c, j;
  float float_tmp;
  int int_tmp;
  // declare the return object
  VISMatrix matrix;
  // get the number rows and cols
  boolean good=true;
  if (VPF::set(r, pf[name][0]) != VPF::VALID)
    good=false;
  if (VPF::set(c, pf[name][1]) != VPF::VALID)
    good=false;
  if (!good)
    {
      cout << "ImplicitModel::bad matrix " << name << endl;
      return VISMatrix();
    }
  // create the return object
  matrix = VISMatrix(r,c);
  // read the data
  for (j = 0; j < r*c; j++)
    if (VPF::set(float_tmp, pf[name][j+2]) == VPF::VALID)
      matrix.poke(j/c, j%c) = float_tmp;
    else 
      if (VPF::set(int_tmp, pf[name][j+2]) == VPF::VALID)
	matrix.poke(j/c, j%c) = (float)int_tmp;
      else
	good=false;

  if (!good)
    {
      cout << "ImplicitModel::bad matrix " << name << endl;
      return VISMatrix();
    }
  return(matrix);
}

// Given an image in float array this function returns a VISImage equivalent object
VISImage<float> getVISImage(float *I, int w, int h) {
	VISImage<float> retI(w,h);
	
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			retI.poke(i,j) = I[j*w+i];
	
	return retI;
}

// Given a VISImage object this function returns a float array equivalent
float* getFloatArray(VISImage<float> I) {
	int w = I.width();
	int h = I.height();
	float* retI = (float*)malloc(sizeof(float)*w*h);
	
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			retI[j*w+i] = I.peek(i,j);
	
	return retI;
}

// Writes the image in the given filenameS
void writeImage(VISImage<float> I, char* filename) {
	VISImageFile im_file;
	
	string of(filename);
	if(endswith(of, ".tif") || endswith(of,".tiff") || endswith(of, ".TIF") || endswith(of, ".TIFF"))
		im_file.write_tiff(VISImage<byte>(I), filename);
	else im_file.write(I, filename);
}
