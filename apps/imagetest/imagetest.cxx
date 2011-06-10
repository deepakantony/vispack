#include <image.h>
#include <imageRGBA.h>
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
  VISImageRGBA image_color;
  VISImage<byte> image_byte;

  // declare an object to read files and produce im objects
  VISImageFile im_file;

  float the_min, the_max;
  int h, w;

  if ((argc > 1)&&(im = im_file.read(argv[1])).isValid())
    {
      // construct an image of type float (this will cast each pixel)
      image = VISImage<float>(im);
      image_color = VISImageRGBA(im);
      image_byte = VISImage<byte>(im);
      image_name = argv[1];
    }
  else
    {
      // make an image of type float of specific size
      image = VISImage<float>(IMAGE_SIZE, IMAGE_SIZE);
      image.evaluate(sinusoid);
      image_name = test_name;
    }

  im_file.write(image_color, "imagetest_color.tif");
  im_file.write_tiff(image_byte, "imagetest_byte.tif");
    
  the_min = image.min(); the_max = image.max();
  cout << "The min and max of " << image_name << " is " << the_min
       << " and " << the_max << endl;

  // get the height and width
  h = image.height(); w = image.width();

  // create a new image of same size
  VISImage<float> image2(w, h);

  // loop through the image and get pixels and square them and put them in the new image
  int i, j;
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	// take the square of the image
	image2.poke(i, j) = image.peek(i, j)*image.peek(i, j);
	// could do this outside this loop with image2 = image.power(2);
      }

  // get the height and width of the input image
  // create a little bit bigger image
  // notice, you don't need template parameter since float is default
  image_out = VISImage<>(w +10, h+10);
  // set the image to all of the same value
  image_out = (the_min + the_max)/2.0f;
  // place the old image into the new one
  image_out.putROI(image, 5, 5);

  // save the new image
  im_file.write_jpeg(image_out, "imagetest_out1.jpg", 100);

  // blur with a gaussian kernel of 2.0  
  // uses zero as the boundary conditions
  image_out = image.gauss(4.0f);
  // to overwrite an old image, you need an "!" before the filename in fits
  im_file.write(image_out, "!imagetest_out2.fits");

  cout << "done." << endl;
  
}

