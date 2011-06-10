#include "image/image.h"
#include "image/imagefile.h"
#include "util/mathutil.h"
#include <stdio.h>
#include <stdlib.h>

// This is a routine from numerical recipes found in the file fourn.c
void fourn(float data[], unsigned long nn[], int ndim, int isign);

// a circularly symmetric butterworth LP filter.
VISImage<float> butterWorthLP(int w, int h, float D0, float n)
{
  VISImage<float> ret(w, h);
  int i, j;
  float dist;
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	dist = power((i - w/2), 2) + power((j - h/2), 2);
	ret.poke(i, j) = 1.0f/(1.0 + exp(n*log(dist/(D0*D0))));
      }
  return(ret);
}


// Does a fourier transform of an image.  Returns two images, as a
// VISArray.  The first one is the real part and the second one the
// imaginary part. 
VISArray< VISImage<float> > fourier(const VISImage<float> &im)
{
  int h_old = im.height(), w_old = im.width();
  int h = pow(2, (int)ceil(log(h_old)/log(2.0))), w = pow(2, (int)ceil(log(w_old)/log(2.0)));
  cout << "w and h" << w_old << " " << h_old << " " << endl;
  cout << "w and h" << w << " " << h << " " << endl;

  VISImage<float> real(w, h), imag(w, h);
  VISImage<float> im_new(w, h);
  im_new = 0.0f;
  int diff_x = w - w_old, diff_y = h - h_old; 
  im_new.putROI(im, diff_x/2, diff_y/2);
  int i, j;

  float* buf_comp = new float[h*w*2];  

  unsigned long dim[2];  
  dim[1] = w; dim[0] = h;

  float *im_ptr = buf_comp;
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	*im_ptr++ = im_new.peek(i, j);
	*im_ptr++ = 0.0;
      }

  fourn(buf_comp - 1, dim - 1, 2, 1);

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	real.poke(i, j) = buf_comp[2*(i*h + j)];
	imag.poke(i, j) = buf_comp[2*(i*h + j) + 1];
      }

  VISArray< VISImage<float> > ret;

  ret.poke(0) = real;
  ret.poke(1) = imag;

  return(ret);
}

// Does the inverse Fourier transform, and returns a real image (the
// real part).    Input is an array of images.  The first image is the
// real part and the second one is the imaginary part. 
// 
VISImage<float> fourierInv(const VISArray< VISImage<float> > &imf)
{
  VISImage<float> real = imf.peek(0), imag = imf.peek(1);
  int w = imag.width(), h = imag.height();
  int i, j;

  float* buf_comp = new float[h*w*2];  

  unsigned long dim[2];  
  dim[1] = w; dim[0] = h;

  float* im_ptr = buf_comp;
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	*im_ptr++ = real.peek(i, j);
	*im_ptr++ = imag.peek(i, j);
      }

  fourn(buf_comp - 1, dim - 1, 2, -1);
  VISImage<float> ret(w, h);

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	ret.poke(i, j) = buf_comp[2*(i*h + j)];
      }

  return(ret/(float)(w*h));

}


// This simply does the complex multiplication of F and G and does the inverse
VISImage<float> filterF(const VISArray< VISImage<float> > &F, 
		       const VISArray< VISImage<float> > &G)
{
  VISArray< VISImage<float> > H;
  (F.peek(0)).print();
  (F.peek(1)).print();
  H.poke(0) = F.peek(0)*G.peek(0) - F.peek(1)*G.peek(1);
  H.poke(1) = F.peek(0)*G.peek(1) + F.peek(1)*G.peek(0);
  return(fourierInv(H));
}

// This reoganized the fourier data with the origin in the middle of
// the image.  
VISImage<float> retile(const VISImage<float> &im)
{
  int w = im.width(), h = im.height();
  VISImage<float> ret(w, h);
  int x, y;
  int i, j;

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	ret.poke(i, j) =  im.peek((i + w/2)%w, (j + h/2)%h);
      }
  return(ret);
}


main(int argc, char **argv)
{
    VISImageFile im_file;
    VISImage<float> im;

    im = VISImage<float>(im_file.read(argv[1]));
    VISArray< VISImage<float> > im_f;
    VISImage<float> out_f;
    im_f = fourier(im);
    out_f = fourierInv(im_f);
    int w, h;
    w = out_f.width(); h = out_f.height();

    cout << "orig min " << im.min() << " and max " << im.max() << endl;
    cout << "f min " << (out_f.min())/(float)(w*h) << " and max " << (out_f.max())/(float)(w*h) << endl;
    
    im_file.write(out_f, "test.fts");

    w = (im_f.peek(0)).width(); h = (im_f.peek(0)).height();
    VISImage<float> filter(w, h), imag(w, h);
    VISArray< VISImage<float> > filter_f;

    imag = 0.0f;
    filter_f.poke(1) = imag;

    filter_f.poke(0) = retile((filter = butterWorthLP(w, h, 30.0, 1.0f)));
    im_file.write(filterF(filter_f, im_f), "butterworth30-1.fit");
    im_file.write(filter, "butterworthfilter30-1.fit");
    filter_f.poke(0) = retile((filter = butterWorthLP(w, h, 30.0, 2.0f)));
    im_file.write(filterF(filter_f, im_f), "butterworth30-2.fit");
    im_file.write(filter, "butterworthfilter30-2.fit");
    filter_f.poke(0) = retile((filter = butterWorthLP(w, h, 30.0, 4.0f)));
    im_file.write(filterF(filter_f, im_f), "butterworth30-4.fit");
    im_file.write(filter, "butterworthfilter30-4.fit");

    cout << "done" << endl;
}
