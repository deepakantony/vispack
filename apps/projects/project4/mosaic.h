#ifndef _MOSAIC_H
#define _MOSAIC_H

#include <imagefile.h>
#include <matrix.h>
#include <image.h>
//#include <cmath>
#include <string>
#include "transform.h"

using namespace std;

class ImageMosaic {
	public:
		ImageMosaic(VISArray< VISImage< float > > images, 
					int num_images,
					VISArray< VISArray< VISArray<int> > > translations,
					VISMatrix translation_matrix,
					string output_filename);
	
	private:
		void initTransformations();
		void computeMosaic();
		CascadeWarp* getTransformationSequence(int imageNum);
		void initFinalImageSize();
		void dijkstras();
		VISArray<float> getPerspectiveParameters(int, int);
		void correctContrast();
		
	private: 
		VISArray<int> shortestPath;
		// the only transformations required to create the image mosaic
		VISArray< CascadeWarp* > T; 
		// the inverse transformations
		VISArray< CascadeWarp* > invT;
		VISMatrix translation_matrix;
		// a 2D array of all the transformations
		VISArray< VISArray< WarpType< Perspective > * > > transformations;
		// a 2D array of the translations - this will be used to create the transformations
		VISArray< VISArray< VISArray<int> > > translations;
		int num_images;
		int target_image;
		VISArray< VISImage< float > > images;
		VISImage< float > finalImage;
		int minx, maxx, miny, maxy;
		int offsetx, offsety;
		string output_filename;
};

#endif /* MOSAIC_H */

