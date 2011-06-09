#include <image.h>
#include <imagefile.h>


#define IMAGE_SIZE 200

// this is a stand along scalar function that we can plug into an image
// to make a bumpy pattern
float sinusoid(unsigned u, unsigned v)
{
  return(sin((2.0*M_PI/IMAGE_SIZE)*(1.0f/10.0f)*(float)u)
	 + sin((2.0*M_PI/IMAGE_SIZE)*(1.0f/20.0f)*(float)v));
}


main(int argc, char** argv)
{

  char *test_name = "testimage";
  char* image_name;

  // a typeless base class
  VISIm im;

  // declare an image of type float --- default size is zero
  VISImage<float> image, image_out;

  // declare an object to read files and produce im objects
  VISImageFile im_file;

  float the_min, the_max;
  int h, w;

  if ((argc > 1)&&(im = im_file.read(argv[1])).isValid())
    {
      // construct an image of type float (this will cast each pixel)
      image = VISImage<float>(im);
      image_name = argv[1];
    }
  else
    {
      // make an image of type float of specific size
      image = VISImage<float>(IMAGE_SIZE, IMAGE_SIZE);
      image.evaluate(sinusoid);
      image_name = test_name;
    }
    
  the_min = image.min(); the_max = image.max();
  cout << "The min and max of " << image_name << " is " << the_min
       << " and " << the_max << endl;

  // get the height and width of the input image
  h = image.height(); w = image.width();
  // create a little bit bigger image
  // notice, you don't need template parameter since float is default
  image_out = VISImage<>(w +10, h+10);
  // set the image to all of the same value
  image_out = (the_min + the_max)/2.0f;
  // place the old image into the new one
  image_out.putROI(image, 5, 5);

  // save the new image
  im_file.write(image_out, "imagetest_out1.fits");

  // blur with a gaussian kernel of 2.0  
  // uses zero as the boundary conditions
  image_out = image.gauss(4.0f);
  im_file.write(image_out, "imagetest_out2.fits");

  // blur with a gaussian kernel of 3.0  
  // uses adiabatic boundary conditions
  image_out = image.gaussDiffuse(3.0f);
  im_file.write(image_out, "imagetest_out3.fits");

  // find the gradient magnitude (edges) of a slightly blurred version
  image_out = image.gaussDiffuse(0.5);
  image_out = ((image.dx()).power(2) + (image.dy()).power(2)).sqrt();
  // set the border to zero, otherwise the boundaries allways look like a strong edge
  image_out = image_out.setBorder(0.0f, 1);
  im_file.write(image_out, "imagetest_out4.fits");

  cout << "done." << endl;
  
}

