#include <image.h>
#include <imagefile.h>



VISImage<float> waveletX(const VISImage<float> image, const VISImage<float> kernel)
{
  int w = image.width(), h = image.height();
  int step = kernel.width();
  if (kernel.height() != 1)
    cout << "warning: x kernel has height" << endl;
  int w_new = w/step;
  VISImage<float> ret(w_new, h);
  float total;
  int i, j, k;
  for (j = 0; j < h; j++)
    for (i = 0; i< w_new; i++)
      {
	total = 0.0f;
	for (k = 0; k < step; k++)
	  total += image.peek(step*i + k, j)*kernel.peek(k, 0);
	ret.poke(i, j) = total;
      }
  cout << "ret min " << ret.min() << " and max " << ret.max() << endl;
  ret.print();
  return(ret);
}

VISImage<float> waveletY(const VISImage<float> image, const VISImage<float> kernel)
{
  int w = image.width(), h = image.height();
  int step = kernel.height();
  if (kernel.width() != 1)
    cout << "warning: x kernel has height" << endl;
  int h_new = h/step;
  VISImage<float> ret(w, h_new);
  float total;
  int i, j, k;
  for (j = 0; j < h_new; j++)
    for (i = 0; i< w; i++)
      {
	total = 0.0f;
	for (k = 0; k < step; k++)
	  total += image.peek(i, step*j + k)*kernel.peek(0, k);
	ret.poke(i, j) = total;
      }
  cout << "ret min " << ret.min() << " and max " << ret.max() << endl;
  ret.print();
  return(ret);
}


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

  cout << "input min " << image.min() << " and max " << image.max() << endl;
  VISImage<float> harrX_lo(2, 1), harrX_hi(2, 1),  harrY_lo(1, 2), harrY_hi(1, 2);
  harrX_lo = 0.5f;
  harrX_hi.poke(0, 0) = 1.0;
  harrX_hi.poke(1, 0) = -1.0;
  harrY_lo = 0.5f;
  harrY_hi.poke(0, 0) = 1.0;
  harrY_hi.poke(0, 1) = -1.0;
  
  int levels = 1;
  VISImage<float> LL, LH, HL, HH;
  LL = image;
  image_out = image.createToSize();
  int positionX, positionY;
  for (int k = 0; k < levels; k++)
    {
      HH = waveletX(waveletY(LL, harrY_hi), harrX_hi);
      HL = waveletX(waveletY(LL, harrY_lo), harrX_hi);
      LH = waveletX(waveletY(LL, harrY_hi), harrX_lo);
      LL = waveletX(waveletY(LL, harrY_lo), harrX_lo);
      positionX = LL.width(); positionY = LL.height();
      image_out.putROI(HH, positionX, positionY);
      image_out.putROI(HL, positionX, 0);
      image_out.putROI(LH, 0, positionY);
    }
  image_out.putROI(LL, 0, 0);
  im_file.write(image_out, "wavelets.fts");
}

