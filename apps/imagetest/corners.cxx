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


VISImage<float> computeCorners(const VISImage<float> &im)
{
  VISImage<float> im_dx, im_dy, nx, ny; 
  VISImage<float> kernel_dx(3,1), kernel_dy(1,3);
  VISImage<float> grad_mag, curve;
  kernel_dx.poke(0, 0) = 0.0f;
  kernel_dx.poke(1, 0) = 0.5f;
  kernel_dx.poke(2, 0) = 0.5f;
  //
  kernel_dy.poke(0, 0) = 0.0f;
  kernel_dy.poke(0, 1) = 0.5f;
  kernel_dy.poke(0, 2) = 0.5f;
  //
  im_dx = im.dx();
  im_dy = im.dy();

  grad_mag = (im_dx.power(2) + im_dy.power(2)).sqrt();

  curve = (im.dx(2)*im_dy.power(2) + 
    im.dy(2)*im_dx.power(2)
	   - im_dx.dy()*im_dx*im_dy).abs()/(grad_mag.power(2) + (float)1.0e-2);
  curve = curve.setBorder(0.0f, 2);
  return(grad_mag*curve.power(2));
  //  grad_mag = (im_dx.power(2) + im_dy.power(2)).sqrt();
  grad_mag = (im_dx.power(2) + im_dy.power(2));
  grad_mag = grad_mag.setBorder(0.0f, 1);

  cout << "done 0" << endl;

  //  im_dx = (im_dx).convolve(kernel_dx);
  //  im_dy = (im_dy).convolve(kernel_dy);
  
  int w = im.width(), h = im.height();
  int i, j;
  nx = im.dxHalfForward();
  //  for (i = 0; i < h; i++)
  //    nx.poke(w-1, i) = 0.0f;
  ny = im.dyHalfForward();
  //  for (i = 0; i < w; i++)
  //    ny.poke(i, h-1) = 0.0f;
  float epsilon = 1.0e-2;
  nx /= (im_dx.power(2) + nx.power(2) + epsilon).sqrt();
  ny /= (im_dy.power(2) + ny.power(2) + epsilon).sqrt();
  return((nx.dxHalfBack() + ny.dyHalfBack()));
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

  image = image.gaussDiffuse(3.0);
  im_file.write(computeCorners(image), "corners.fit");
  cout << "done" << endl;
}

