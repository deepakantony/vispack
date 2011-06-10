#include <imagefile.h>
#include "edgedetector.h"
#include "util.h"


int main(int argc, char* argv[]) {
	if(argc < 4) {
		cout << "Usage: <program> <image> <aniso-conductance> <aniso-iterations> <output_filename>" << endl;
		return -1;
	}
	
	VISIm im;
	VISImageFile imgFile;
	VISImage<float> I;
	if((im = imgFile.read(argv[1])).isValid()) {
		I = VISImage<float>(im);
	}
	
	// Pass it through a median filter
	float its;
	float conductance;
	sscanf(argv[2], "%f", &conductance);
	sscanf(argv[3], "%f", &its);
	VISImage<float> aniso = I;
	for(int i = 0; i < its; i++)
		aniso = aniso.anisoDiffuse(conductance);
	
	if(argc > 4)
		writeImage(aniso, argv[4]);
	else writeImage(aniso, "edge.tif");
}
