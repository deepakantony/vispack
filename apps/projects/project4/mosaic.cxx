#include "mosaic.h"


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

ImageMosaic::ImageMosaic (VISArray< VISImage< float > > images, 
						  int num_images,
						  VISArray< VISArray< VISArray<int> > > translations, 
						  VISMatrix translation_matrix,
						  string output_filename) 
{
	this->images = images;
	this->translations = translations;
	this->translation_matrix = translation_matrix;
	this->output_filename = output_filename;
	this->num_images = num_images;
	this->target_image = 0;
								  
	// Initialize the transformations
	initTransformations();
	
	// Contrast correction
	correctContrast();
	
	// Find the offset and resolution of the final image
	initFinalImageSize();
	
	// Find the image mosaic
	computeMosaic();
}

void ImageMosaic::initTransformations() {
	cout << "Initializing the transformations" << endl;
	
	// Loop through the correspondence matrix and set up the perspective projections.
	int size = images.n();
	for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++)
			if(translation_matrix.peek(i,j)) {
				cout << "Creating perspective transformation from image " << i << " to image " << j << endl;
				
				// Transformation from i to j
				VISArray<float> pij(getPerspectiveParameters(i, j));
				pij.poke(8) = i;
				pij.poke(9) = j;
				WarpType<Perspective> *persIJ = new WarpType<Perspective>();
				persIJ->setParams(pij);
				transformations.poke(i).poke(j) = persIJ;
				
				cout << "Creating perspective transformation from image " << j << " to image " << i << endl;
				
				// Transformation from j to i
				VISArray<float> pji(pij);
				pji.poke(8) = j;
				pji.poke(9) = -1; // indicates inverse transformation
				WarpType<Perspective> *persJI = new WarpType<Perspective>();
				persJI->setParams(pji);
				transformations.poke(j).poke(i) = persJI;
			}
	
	cout << "Done setting up transformations from the initial values" << endl;
	
	// Add its own transpose as we have set up the transformations for i->j and j->i transformations
	translation_matrix += translation_matrix.t();
	
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
VISArray<float> ImageMosaic::getPerspectiveParameters(int i, int j) {
	VISArray<int> pos = translations.poke(i).poke(j);
	
	VISArray<float> result;
	
	// diagonals
	result.poke(0) = 1.f;
	result.poke(4) = 1.f;
	result.poke(8) = 1.f;
	
	// translation component
	result.poke(2) = pos.peek(0);
	result.poke(5) = pos.peek(1);
	
	// others are 0!
	result.poke(1) = 0.f;
	result.poke(3) = 0.f;
	result.poke(6) = 0.f;
	result.poke(7) = 0.f;
	
	return result;
	
/*	cout << "Setting up the system Ax = b." << endl;

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
	
	return result;*/
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
		im_file.write_tiff(VISImage<byte>(finalImage), output_filename.c_str());
	else im_file.write(finalImage, output_filename.c_str());
	
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
			if(Q.peek(i) && translation_matrix.peek(u, i)) {
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
			if(translation_matrix.poke(i,j) && (!correctedImage[i] || !correctedImage[j])) {
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
