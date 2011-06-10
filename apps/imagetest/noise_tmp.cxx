#include <image.h>
#include <imagefile.h>


#define IMAGE_SIZE 200

// this is a stand along scalar function that we can plug into an image
// to make a bumpy pattern

main(int argc, char** argv)
{

  VISImage<float> image;
  VISIm im;
  VISImageFile im_file;
  float the_min, the_max;

  if ((argc > 1)&&(im = im_file.read(argv[1])).isValid())
    {
      // construct an image of type float (this will cast each pixel)
      image = VISImage<float>(im);
    }
    
  the_min = image.min(); the_max = image.max();
  cout << "The min and max is " << the_min
       << " and " << the_max << endl;

  //  image = image.noiseShot(0.2, 
  //			  the_min - 0.1*(the_max - the_min), 
  //			  the_max + 0.1*(the_max - the_min));

  image += 0.1f*(the_max - the_min)*noiseGauss(image.width(), image.height(), 1.0f);

  im_file.write(image, "noise_tmp.tif");
  
}

