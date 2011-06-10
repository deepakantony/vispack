#include <image.h>
#include <imagefile.h>

main(int argc, char** argv)
{

  char filename[80];

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
    }
  else
    {
      cout << "bad input image " << endl;
    }
  im_file.write((image.dx()).setBorder(0.0, 2), "image_dx.fits");
  im_file.write((image.dy()).setBorder(0.0, 2), "image_dy.fits");
  im_file.write(((image.dx()).abs() + (image.dy()).abs()).setBorder(0.0, 2), "image_grad_L1.fits");
  im_file.write((((image.dx()).power(2) + (image.dy()).power(2)).sqrt()).setBorder(0.0, 2), "image_grad_mag.fits");
  im_file.write(image.cannyEdges(10.0), "image_edges.fits");

  image = image.gauss(2.0f);
  im_file.write((image.dx()).setBorder(0.0, 4), "image_2.0_dx.fits");
  im_file.write((image.dy()).setBorder(0.0, 4), "image_2.0_dy.fits");
  im_file.write(((image.dx()).abs() + (image.dy()).abs()).setBorder(0.0, 4), "image_2.0_grad_L1.fits");
  im_file.write((((image.dx()).power(2) + (image.dy()).power(2)).sqrt()).setBorder(0.0, 4), "image_2.0_grad_mag.fits");
  im_file.write(image.cannyEdges(5.0), "image_2.0_edges.fits");

  // find the gradient magnitude (edges) of a slightly blurred version
  //  image_out = image.gaussDiffuse(0.5);

  cout << "done." << endl;
  
}

