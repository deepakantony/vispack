#include <array.h>
#include <param.h>
#include <imagefile.h>
#include <matrix.h>
#include <image.h>
#include <cmath>
#include <string>
#include "transform.h"

using namespace std;

template<class T>
void displayArray(VISArray<T> array) {
	for(int i = 0; i < array.n(); i++)
		cout << array[i] << " ";
	cout << endl;
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

//
// A helper function which reads, from a ParameterFile object, a
// matrix of the form (name rows cols ....data in row major order...)
// Returns an invalid matrix if it fails (prints something)
//
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

class ImageMosaic {
	public:
		ImageMosaic(char *paramfile);
	
	private:
		void initParams(char *paramfile);
		void initTransformations();
		void computeMosaic();
		CascadeWarp* getTransformationSequence(int imageNum);
		void initFinalImageSize();
		void dijkstras();
		VISArray<float> getPerspectiveParameters(VISMatrix, VISMatrix);
		void correctContrast();
		
	private: 
		VISArray<int> shortestPath;
		// the only transformations required to create the image mosaic
		VISArray< CascadeWarp* > T; 
		// the inverse transformations
		VISArray< CascadeWarp* > invT;
		VISMatrix correspondence_matrix;
		// create a 2D array that contains all of the control points
		VISArray< VISArray< VISMatrix > >  correspondences;
		// a 2D array of all the transformations
		VISArray< VISArray< WarpType< Perspective > * > > transformations;
		int num_images;
		int target_image;
		VISArray< VISImage< float > > images;
		VISImage< float > finalImage;
		int minx, maxx, miny, maxy;
		int offsetx, offsety;
		char output_filename[80];
};

ImageMosaic::ImageMosaic (char *paramfile) {
	// Initialize the correspondence matrices and images
	initParams(paramfile);
	
	// Initialize the transformations
	initTransformations();
	
	// Contrast correction
	correctContrast();
	
	// Find the offset and resolution of the final image
	initFinalImageSize();
	
	// Find the image mosaic
	computeMosaic();
}
	
void ImageMosaic::initParams(char *paramfile) {
	// create parameter file from the given argument
	VPF::ParameterFile pf(paramfile);	
	
	// read in the matrix of correspondences for which we will look  
	correspondence_matrix = readMatrix(pf, "CORRESPONDENCE_MATRIX");

	// print it out just for debugging
	cout << "correspondence matrix " << endl << correspondence_matrix << endl;

	// make sure it's valid and square
	if ((!correspondence_matrix.isValid())
		||(!(correspondence_matrix.r() == correspondence_matrix.c()))) {
		cout << "bad correspondence matrix in param file: " << paramfile << endl;
		exit(-1);
	}

	// the number of images is the size of this matrix  (one row/col per
	// image)
	num_images = correspondence_matrix.r();

	int i, j, k;
	char keyword[80];
	VISMatrix pointsA, pointsB;


	// put empty matrices along the diagonal
	for (i = 0; i < num_images; i++)
	(correspondences.poke(i)).poke(i) = VISMatrix();

	// for all entries in the upper half
	for (i = 0; i < num_images; i++)
		for (j = i+1; j < num_images; j++)
			if (correspondence_matrix.peek(i, j)) {
				// look for the two matching correpondence matrices
				sprintf(keyword, "%s_%d_%d", "CORRESPONDENCES", i, j);
				pointsA = readMatrix(pf, keyword);
				sprintf(keyword, "%s_%d_%d", "CORRESPONDENCES", j, i);
				pointsB = readMatrix(pf, keyword);
				// make sure they are the same size
				if (!pointsB.compareSize(pointsA))
				{
				cout << "bad control matrices " << i << " " << j << " in param file: " << paramfile << endl;
				exit(-1);
				}
				// transpose so that individual points are columns
				// and put them in the 2D array of matrices
				(correspondences.poke(i)).poke(j) = pointsA.t();
				(correspondences.poke(j)).poke(i) = pointsB.t();
			}
			else {
				// if those correspondences are not expectd, put
				// in invalid matrices
				(correspondences.poke(i)).poke(j) = VISMatrix();
				(correspondences.poke(j)).poke(i) = VISMatrix();
			}

	// get the target image number  
	if ((VPF::set(target_image, pf["TARGET_IMAGE"][0]) != VPF::VALID)
		||(target_image >= num_images)) {
		cout << "bad target image in param file: " << paramfile << endl;
		exit(-1);
	}

	// print out the correspondeces --- for debugging
	for (i = 0; i < num_images; i++)
		for (j = i+1; j < num_images; j++)
			if (correspondence_matrix.peek(i, j)) {
				cout << "The correspondences between image " << i 
				<< " and image " << j << " are: "  << endl;
				pointsA = correspondences[i][j];
				pointsB = correspondences[j][i];
				for (k = 0; k < pointsA.c(); k++) {
					cout << "c" << k << "=(" << pointsA.peek(0, k) << "," << pointsA.peek(1, k) << ")  : c'" << 
					k << "=(" << pointsB.peek(0, k) << "," << pointsB.peek(1, k) << ") " << endl;
				}
				cout << endl;
			}

	VISImageFile im_file;
	// read in the image and put them in an array
	char filename[80];
	for (i = 0; i < num_images; i++) {
		sprintf(keyword, "%s_%d", "IMAGE_FILE", i);    
		if (VPF::set(filename, pf[keyword][0]) != VPF::VALID) {
			cout << "no image file " << keyword << " in param file: " << paramfile << endl;
			exit(-1);
		}
		images.poke(i) = VISImage<float>(im_file.read(filename));
		if (!(images.peek(i)).isValid()) {
			cout << "bad image file " << filename << " in param file: " << paramfile << endl;
			exit(-1);
		}
	}
	
	// get the output image name
	if (VPF::set(output_filename, pf["OUTPUT_FILE_NAME"][0]) != VPF::VALID) {
		cout << "no output filename specified: " << paramfile << endl;
		exit(-1);
	}
}

void ImageMosaic::initTransformations() {
	cout << "Initializing the transformations" << endl;
	
	// Loop through the correspondence matrix and set up the perspective projections.
	int size = images.n();
	for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++)
			if(correspondence_matrix.peek(i,j)) {
				cout << "Creating perspective transformation from image " << i << " to image " << j << endl;
				
				// Transformation from i to j
				VISArray<float> pij(getPerspectiveParameters(correspondences[i][j].t(), correspondences[j][i].t()));
				pij.poke(8) = i;
				pij.poke(9) = j;
				WarpType<Perspective> *persIJ = new WarpType<Perspective>();
				persIJ->setParams(pij);
				transformations.poke(i).poke(j) = persIJ;
				
				cout << "Creating perspective transformation from image " << j << " to image " << i << endl;
				
				// Transformation from j to i
				VISArray<float> pji(pij);//getPerspectiveParameters(correspondences[j][i].t(), correspondences[i][j].t()));
				pji.poke(8) = j;
				pji.poke(9) = -1; // indicates inverse transformation
				WarpType<Perspective> *persJI = new WarpType<Perspective>();
				persJI->setParams(pji);
				transformations.poke(j).poke(i) = persJI;
			}
	
	cout << "Done setting up transformations from the initial values" << endl;
	
	// Add its own transpose as we have set up the transformations for i->j and j->i transformations
	correspondence_matrix += correspondence_matrix.t();
	
	cout << "Calling dijkstras" << endl;
	// Identify paths between target and every other image.
	// using dijkstras algorithm
	dijkstras();
	
	cout << "Done with executing dijkstras algorithm" << endl;
	cout << "Shortest paths array: ";
	for(int i = 0; i < shortestPath.n(); i++)
		cout << shortestPath.peek(i) << " ";
	cout << endl;
	
	// loop through the paths and build the transformations for all images
	// create inverse transformations along with the forward transformations
	// i.e. build T and invT
	T.sizeTo(size);
	invT.sizeTo(size);
	WarpType< Transform > noTransform;
	T.poke(target_image) = new CascadeWarp(noTransform, noTransform);
	invT.poke(target_image) = new CascadeWarp(*(T.peek(target_image)));
	for(int i = 0; i < shortestPath.n(); i++) {
		if(i == target_image) {
			i++;
			
			if(i >= shortestPath.n()) break;
		}
		int curSP = shortestPath.peek(i);
		
		WarpType<Perspective> firstT(*(transformations.peek(i).peek(curSP)));
		WarpType<Perspective> firstInvT(*(transformations.peek(curSP).peek(i)));
		if(curSP == target_image) {
			T.poke(i) = new CascadeWarp(firstT, noTransform);
			invT.poke(i) = new CascadeWarp(firstInvT, noTransform);
		}
		else {
			int prevSP = curSP;
			curSP = shortestPath.peek(curSP);
			WarpType<Perspective> secondT(*(transformations.peek(prevSP).peek(curSP)));
			WarpType<Perspective> secondInvT(*(transformations.peek(curSP).peek(prevSP)));
			CascadeWarp *prevT = new CascadeWarp(secondT, firstT);
			//CascadeWarp *prevInvT = new CascadeWarp(secondInvT, firstInvT);
			//CascadeWarp *prevT = new CascadeWarp(firstT, secondT);
			CascadeWarp *prevInvT = new CascadeWarp(firstInvT, secondInvT);
			
			prevSP = curSP;
			curSP = shortestPath.peek(curSP);
			while(curSP != -1) {
				secondT = WarpType<Perspective>(*(transformations.peek(prevSP).peek(curSP)));
				secondInvT = WarpType<Perspective>(*(transformations.peek(curSP).peek(prevSP)));
				
				prevT = new CascadeWarp(secondT, *prevT);
				//prevInvT = new CascadeWarp(secondInvT, *prevInvT);
				//prevT = new CascadeWarp(*prevT, secondT);
				prevInvT = new CascadeWarp(*prevInvT, secondInvT);
				
				prevSP = curSP;
				curSP = shortestPath.peek(curSP);
			}
			T.poke(i) = prevT;
			invT.poke(i) = prevInvT;
		}
	}
	
	cout << "Done initializing the required transformations." << endl;
}

// Perspective transformation used correspondences to calculate the perspective transformation.
VISArray<float> ImageMosaic::getPerspectiveParameters(VISMatrix correspondence1, VISMatrix correspondence2) {
	cout << "Setting up the system Ax = b." << endl;

	VISMatrix c2 = correspondence1;
	VISMatrix c1 = correspondence2;
	
	VISVector c1x = c1.vec(0);
	VISVector c1y = c1.vec(1);
	VISVector c2x = c2.vec(0);
	VISVector c2y = c2.vec(1);
	
	VISVector b = c1x;
	b = -b.concatRow(c1y);
	
	cout << "b - done" << endl;
	
	cout << "Building A" << endl;
	VISMatrix temp1(c2.r(), c2.c()+1);
	temp1 = 0;
	VISMatrix temp2(c2.r(), 1);
	temp2 = 1;
	
	c2 = c2.concat(temp2);
	
	VISMatrix A = -c2;
	A = A.concatRow(temp1);

	VISMatrix temp11(c2.r(), 1);
	temp11 = 0;
	temp11 = temp11.concatRow(-c2x);
	VISMatrix temp22(c2.r(), 1);
	temp22 = 0;
	temp22 = temp22.concatRow(-c2y);
	VISMatrix temp33(c2.r(), 1);
	temp33 = 0;
	temp33 = temp33.concatRow(-temp2);
	A = A.concat(temp11.concat(temp22.concat(temp33)));
	
	VISVector temp3 = mul(c1x, c2x);
	temp3 = temp3.concatRow(mul(c1y,c2x));
	A = A.concat(temp3);
	
	VISVector temp4 = mul(c1x, c2y);
	temp4 = temp4.concatRow(mul(c1y, c2y));
	A = A.concat(temp4);
	
	cout << "A - done" << endl;
	
	// Ax = b, solve this linear system
	VISVector x = solvLinearEqn(A, b);
	
	cout << "Linear system solved." << endl;
	
	VISArray<float> result;
	for(int i = 0; i < x.n(); i++)
		result.poke(i) = x.peek(i);
	
	return result;
}

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

// Finds the final images size
// This function transforms the corners of the image to identify the image sizes.
void ImageMosaic::initFinalImageSize() {
	cout << "Finding the image size...." << endl;
	
	int minx = +99999;
	int miny = +99999;
	int maxx = -99999;
	int maxy = -99999;
	
	for(int i = 0; i < images.n(); i++) {
		// Coordinate 0,0
		VISVector c1(2);
		c1.poke(0) = 0;
		c1.poke(1) = 0;
		updateMinsMaxs(T.peek(i)->transform(c1), minx, miny, maxx, maxy);
		
		// Coordinate width, 0
		VISVector c2(2);
		c2.poke(0) = images.peek(i).width()-1;
		c2.poke(1) = 0;
		updateMinsMaxs(T.peek(i)->transform(c2), minx, miny, maxx, maxy);
			
		// Coordinate 0, height
		VISVector c3(2);
		c3.poke(0) = 0;
		c3.poke(1) = images.peek(i).height()-1;
		updateMinsMaxs(T.peek(i)->transform(c3), minx, miny, maxx, maxy);
			
		// Coordinate width, height
		VISVector c4(2);
		c4.poke(0) = images.peek(i).width()-1;
		c4.poke(1) = images.peek(i).height()-1;
		updateMinsMaxs(T.peek(i)->transform(c4), minx, miny, maxx, maxy);
	}
	
	cout << "MinX: " << minx << endl;
	cout << "MinY: " << miny << endl;
	cout << "MaxX: " << maxx << endl;
	cout << "MaxY: " << maxy << endl;
	
	this->minx = minx;
	this->miny = miny;
	this->maxx = maxx;
	this->maxy = maxy;
	this->offsetx = -minx;
	this->offsety = -miny;
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

// Computes the final mosaic image
// Loops through each point and resamples the data with bilinear interpolation
void ImageMosaic::computeMosaic() {
	int w = maxx + offsetx;
	int h = maxy + offsety;
	
	finalImage = VISImage<float>(w,h);
	
	VISArray<VISVector> center;
	for(int k = 0; k < images.n(); k++) {
		VISVector c(2);
		c.poke(0) = ((float)(images[k].width()-1.f))/2.f;
		c.poke(1) = ((float)(images[k].height()-1.f))/2.f;
		center.poke(k) = c;
	}
	
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			float avg = 0.f;
			float totalIntensity = 0.f;
			
			// Loop through the images and find out which images are affecting the mosaic
			for(int k = 0; k < images.n(); k++) {
				VISVector x = VISVector((float)i, (float)j) - VISVector((float)offsetx, (float)offsety);
				
				VISVector x_prime = invT[k]->transform(x);
				
				if(images[k].checkBounds((float)x_prime.peek(0), (float)x_prime.peek(1))) {
					float weight = weightFunction(x_prime, center[k]);
					avg+=weight;
					totalIntensity += weight*images[k].interp((float)x_prime.peek(0), (float)x_prime.peek(1));
				}
			}
			
			if(avg == 0)
				finalImage.poke(i,j) = 0.f;
			else
				finalImage.poke(i,j) = totalIntensity/(float)avg;
		}
	
	VISImageFile im_file;
	
	string of(output_filename);
	if(endswith(of, ".tif") || endswith(of,".tiff") || endswith(of, ".TIF") || endswith(of, ".TIFF"))
		im_file.write_tiff(VISImage<byte>(finalImage), output_filename);
	else im_file.write(finalImage, output_filename);
	
	cout << "Mosaiced image written to " << output_filename << endl;	
}

bool isEmpty(VISArray<bool> A) {
	for(int i = 0; i < A.n(); i++)
		if(A.peek(i))
			return false;
	return true;
}

void ImageMosaic::dijkstras() {
	int infinity = 9999;
	int size = images.n();
	VISArray<int> dist(size);
	shortestPath = VISArray<int>(size);
	
	cout << "Dijkstras Initialization." << endl;
	// Initializations
	for(int i = 0; i < size; i++) {
		dist.poke(i) = infinity;
		shortestPath.poke(i) = -1;
	}
	dist.poke(target_image) = 0;
	
	VISArray<bool> Q(size);
	for(int i = 0; i < size; i++)
		Q.poke(i) = true;

	cout << "Dijkstras main loop." << endl;
	while(!isEmpty(Q)) {
		int u;
		int min = infinity;
		for(int i = 0; i < size; i++)
			if(Q.peek(i) && dist.peek(i) < min) {
				u = i;
				min = dist.peek(i);
			}
		// u removed from Q
		Q.poke(u) = false;
		
		VISArray<int> neighbors;
		int index = 0;
		for(int i = 0; i < size; i++)
			if(Q.peek(i) && correspondence_matrix.peek(u, i)) {
				neighbors.poke(index) = i;
				index++;
			}
		
		for(int i = 0; i < neighbors.n(); i++) {
			int newdist = dist.peek(u) + 1;
			
			if(newdist < dist.peek(neighbors.peek(i))) {
				dist.poke(neighbors.peek(i)) = newdist;
				shortestPath.poke(neighbors.peek(i)) = u;
			}
		}
	}
}

// This function adjusts the contrasts in the initial image to
void ImageMosaic::correctContrast() {
	cout << "Starting contrast correction" << endl;
	
	VISArray<bool> correctedImage;
	int size = images.n();
	for(int i = 0; i < size; i++) {
		correctedImage.poke(i) = false;
	}

	// Loop through the existing correspondences and find the overlapping region between those images.
	// Find the difference in mean intensity in the overlapping region. Add the difference to the other image.
	for(int i = 0; i < size; i++) {
		for(int j = i+1; j < size; j++) {
			if(correspondence_matrix.poke(i,j) && (!correctedImage[i] || !correctedImage[j])) {
				WarpType<Perspective>* T = transformations[i][j];
				
				float mean1 = 0.f;
				float denom1 = 0.f;
				float mean2 = 0.f;
				float denom2 = 0.f;
				
				// So I need to transform from i'th image to j'th
				int w = images[i].width();
				int h = images[i].height();
				for(int u = 0; u < w; u++)
					for(int v = 0; v < h; v++) {
						VISVector x((float)u, (float)v);
						VISVector x_prime = T->transform(x);
						
						if(images[j].checkBounds((float)x_prime.peek(0), (float)x_prime.peek(1))) {
							mean1 += images[i].peek(u,v);
							denom1 += 1.f;
							mean2 += images[j].peek(x_prime.peek(0), x_prime.peek(1));
							denom2 += 1.f;
						}
					}
				
				mean1 = mean1/denom1;
				mean2 = mean2/denom2;
				
				// Correct the image thats not yet corrected
				if(correctedImage[i]) {
					images.poke(j) += (mean1 - mean2);
					
					correctedImage.poke(j) = true;
				}
				else if(correctedImage[j]) {
					images.poke(i) += (mean2 - mean1);
					
					correctedImage.poke(i) = true;
				}
				else if((mean1-mean2) > 0) {
					images.poke(j) += (mean1 - mean2);
					
					correctedImage.poke(j) = true;
					correctedImage.poke(i) = true;
				}
				else {
					images.poke(i) += (mean2 - mean1);
					
					correctedImage.poke(i) = true;
					correctedImage.poke(j) = true;
				}
			}
		}
	}
	cout << "Ending contrast correction" << endl;
}

//
// Program reads in a parameter file to specify a mosaicing operation
// in accordance with Project 2.  The format of the file is ...
//
//CORRESPONDENCE_MATRIX <int> <int> <floats in row major order> Matrix
//(square) --- rows and cols are image numbers, and upper half of this
//matrix gives correspondences available //IMAGE_FILE_<N> <string>
//Gives image file name for the Nth image
//
// CORRESPONDENCES_<M>_<N> <int> 2 <floats in row major order
// Specifies a Px2 matrix (p is # of correspondences) which give the
// control points in the Mth image that will match up with those in the 
// Nth image.  Note there must be a CORRESPONDENCES_<N>_<M> and it 
// must be the same size
// 
// TARGET_IMAGE <int> specify which target image coordinate system to
// produce the final result
int main(int argc, char** argv)
{
	if (argc < 2) {
		cout << "usage: readparams parameterfile" << endl;
		exit(-1);
	}


	// create an image mosaic object
	ImageMosaic im(argv[1]);
	
}
