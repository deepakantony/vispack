#include <image.h>
#include <imagefile.h>


#define IMAGE_SIZE 200





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

  VISImage<float> dx, dy, dxx, dxy, dyy, grad_mag_sq; 
  
  for (int k = 0; k < 40; k++)
    {
      dx = image.dx(); 
      dy = image.dy(); 
      dxy = dx.dy();
      dxx = image.dx(2);
      dyy = image.dy(2);
      grad_mag_sq = dx.power(2) + dy.power(2) + (float)1.0e-5;
      image += 0.1f*(dxx + dyy - ((dxx*dx.power(2) + dyy*dy.power(2) + 2.0f*dx*dy*dxy))/grad_mag_sq);
    }
  
  im_file.write(image, "tv.fts");
  cout << "done" << endl;
}

