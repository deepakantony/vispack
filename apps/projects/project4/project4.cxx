#include "image/image.h"
#include "image/imagefile.h"
#include "util/mathutil.h"
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>

#include "mosaic.h"

int num_images;

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

// Just does the complem multiplication
VISArray< VISImage<float> > cmplxMult(VISArray< VISImage<float> > F, VISArray< VISImage<float> > G) {
	VISArray< VISImage<float> > H;
 	H.poke(0) = F.peek(0)*G.peek(0) - F.peek(1)*G.peek(1);
	H.poke(1) = F.peek(0)*G.peek(1) + F.peek(1)*G.peek(0);
	return H;
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

// this function calculates the phase correlation
VISImage<float> phaseCorrelation(VISImage<float> f, VISImage<float> g) {
	VISArray< VISImage<float> > F = fourier(f);
	VISArray< VISImage<float> > G = fourier(g);
	
	// Find conjugate
	VISArray< VISImage<float> > F_conj;
	F_conj.poke(0) = F.peek(0);
	F_conj.poke(1) = (F.peek(1) * -1.f);
	
	int w = F.peek(0).width();
	int h = F.peek(0).height();
	
	// lowpass filter
	VISImage<float> lowpass = retile(butterWorthLP(w,h,50.f,1));
	
	// write this image for reporting
	VISImageFile im_file;
	//im_file.write(butterWorthLP(w,h,50.f,1), "bw_50_1.fts");
	
	VISImage<float> real(w,h);
	VISImage<float> imag(w,h);
	
	VISImage<float> realF = F_conj.peek(0);
	VISImage<float> imagF = F_conj.peek(1);
	VISImage<float> realG = G.peek(0);
	VISImage<float> imagG = G.peek(1);
	
	VISImage<float> magF = (sqr(realF)+sqr(imagF)).sqrt();
	VISImage<float> magG = (sqr(realG)+sqr(imagG)).sqrt();
	VISImage<float> mag = magF*magG;
	
	VISArray< VISImage<float> > H;
	H.poke(0) = lowpass;
	H.poke(1) = VISImage<float>(w,h);
	H.poke(1) = 0.f;
	
	// calculate the phase correlation
	VISArray< VISImage<float> > ret = cmplxMult(H, cmplxMult(F_conj, G));
	ret.poke(0) = ret.peek(0)/mag;
	ret.poke(1) = ret.peek(1)/mag;
	
	return fourierInv(ret);
	
	//VISArray< VISImage<float> > H = cmplxMult(F_conj, G);
	
	//real = lowpass*H.peek(0)/mag;
	//imag = lowpass*H.peek(1)/mag;
	
	//real = (lowpass*(F_conj.poke(0)*G.poke(0)-F_conj.poke(1)*G.poke(1)))/((sqr(F_conj.poke(0))+sqr(F_conj).sqrt()*(sqr(G.poke(0))).sqrt());
	//imag = (lowpass*(F_conj.poke(1)*G.poke(0)+F_conj.poke(0)*G.poke(1)))/((sqr(F_conj.poke(1))).sqrt()*(sqr(G.poke(1))).sqrt());
	//real = (lowpass*(F_conj.poke(0)*G.poke(0)))/((sqr(F_conj.poke(0))).sqrt()*(sqr(G.poke(0))).sqrt());
	//imag = (lowpass*(F_conj.poke(1)*G.poke(1)))/((sqr(F_conj.poke(1))).sqrt()*(sqr(G.poke(1))).sqrt());
/*	
	for(int u = 0; u < w; u++) {
		for(int v = 0; v < h; v++) {
			float realF = F_conj.peek(0).peek(u,v);
			float imagF = F_conj.peek(1).peek(u,v);
			float realG = G.peek(0).peek(u,v);
			float imagG = G.peek(1).peek(u,v);
			float magF = sqrt(realF*realF+imagF*imagF);
			float magG = sqrt(realG*realG+imagG*imagG);
			float mag = magF*magG;
			real.poke(u,v) = lowpass(u,v)*(realF*realG-imagF*imagG)/mag;
			imag.poke(u,v) = lowpass(u,v)*(realF*imagG+imagF*realG)/mag;
		}
	}
	
	VISArray< VISImage<float> > pc;
	pc.poke(0) = real;
	pc.poke(1) = imag;
	
	return fourierInv(pc);*/
}

template <class T>
float *histogram(VISImage<T> I, int n, T min, T max) {
	float *H = (float *)malloc(sizeof(float)*n);
	
	for(int i = 0; i < n; i++)
		H[i] = 0.0f;
	
	T mulfact = (n-1)/(max-min);
	int h = I.height();
	int w = I.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			T intensity = I.peek(i,j);
			if( intensity >= min && intensity <= max)
				H[(int)((intensity-min)*mulfact)] += 1.0f;
		}
	
	return H;
}

template<class T>
VISImage<byte> threshold(VISImage<T> I, T min, T max) {
	VISImage<byte> outputImg = I.createToSize();
	outputImg = (T)0.0;
	int h = I.height();
	int w = I.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			if(I.peek(i,j) >= min && I.peek(i,j) <= max)
				outputImg.poke(i,j) = 255;
	
	return outputImg;
}

VISArray< VISImage<float> > pad(VISArray< VISImage<float> > images) {
	int maxw = 0;
	int maxh = 0;
	
	// Find the maximum size and enlarge the image to twice this size
	for(int i = 0; i < num_images; i++) {
		int h_old = images.peek(i).width(), w_old = images.peek(i).height();
  		int h = pow(2, (int)ceil(log(h_old)/log(2.0)))*2, w = pow(2, (int)ceil(log(w_old)/log(2.0)))*2;
		
		if(maxw < w) maxw = w;
		if(maxh < h) maxh = h;
	}
	
	VISArray< VISImage<float> > ret;
	
	// Pad the extra pixels on the right and bottom with average values
	for(int i = 0; i < num_images; i++) {
		VISImage<float> ithIMG(maxw, maxh);
		VISImage<float> curIMG = images.peek(i);
		int curw = curIMG.width();
		int curh = curIMG.height();
		float avg = curIMG.average();
		for(int u = 0; u < maxw; u++)
			for(int v = 0; v < maxh; v++) {
				if(u < curw && v < curh)
					ithIMG.poke(u,v) = curIMG.peek(u,v);
				else
					ithIMG.poke(u,v) = avg;
			}
		
		ret.poke(i) = ithIMG;
	}
	
	cout << "Padding done" << endl;
	
	return ret;
}

VISArray<int> findPeak(VISImage<float> f, VISImage<float> g, float threshold) {
	// get the phase correlation image
	VISImage<float> pc = phaseCorrelation(f,g);
	
	cout << "PC - max: " << pc.max() << ", PC - min: " << pc.min() << endl;
	cout << "PC - average: " << pc.average() << endl;
	
	// Find the peak position
	VISArray<int> ret;
	ret.poke(0) = -1;
	ret.poke(1) = -1;
	
	// temporary variable to store the maximum intensity
	float max = -99999;
	int w = pc.width(), h = pc.height();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			if(pc.peek(i,j) > threshold && pc.peek(i,j) > max) {
				max = pc.peek(i,j);
				ret.poke(0) = i;
				ret.poke(1) = j;
			}
	
	cout << "Position: (" << ret.peek(0) << ", " << ret.poke(1) << ")" << endl;
	
	return ret;
}

void mosaic(VISArray< VISImage<float> > images, string output_filename) {
	// get the threshold value
	float threshold;
	cout << "Enter threshold: ";
	cin >> threshold;
	
	// pad the images
	VISArray< VISImage<float> > newImages = pad(images);
	
	// translation matrix: stores the existing translations
	VISMatrix translation_matrix(num_images, num_images);
	translation_matrix = 0;
	
	// translations: the (X,Y) translation needed from image i to j
	VISArray< VISArray< VISArray<int> > > translations;
	
	// Find the phase correlation peak between each pair
	for(int i = 0; i < (num_images-1); i++)
		for(int j = i+1; j < num_images; j++) {
			// Find the peak for correlation between i and j
			cout << i << " to " << j << " : find peak" << endl;
			VISArray<int> peakPos = findPeak(newImages.peek(i), newImages.peek(j), threshold);
			if(peakPos.peek(0) >= 0) {
				translation_matrix.poke(i,j) = 1;
				
				int w = images.peek(j).width();
				int h = images.peek(j).height();
				
				// padded size
				int pw = newImages.peek(j).width();
				int ph = newImages.peek(j).height();
				
				int x = peakPos.peek(0);
				int y = peakPos.peek(1);
				
				if(x >= w)
					x = x - pw;
				if(y >= h)
					y = y - ph;
				
				// set up translations from i to j
				VISArray<int> pos;
				pos.poke(0) = x; pos.poke(1) = y;
				
				translations.poke(i).poke(j) = pos;
			}
		}
	
	cout << "Translation matrix:\n" << translation_matrix << endl;
	
	// construct the mosaic
	ImageMosaic mosaic(images, num_images, translations, translation_matrix, output_filename);
}

main(int argc, char **argv)
{
	VISImageFile im_file;
	VISArray< VISImage<float> > images;
	
	num_images = argc-1;
	for(int i = 1; i < argc; i++) {
		cout << "Reading " << argv[i] << endl;
		images.poke(i-1) = VISImage<float>(im_file.read(argv[i]));
	}
	
	mosaic(images, "output_mosaic2.fts");
	
	/*images = pad(images);
	VISImage<float> pc = phaseCorrelation(images.peek(0), images.peek(1));
	pc.print();
	cout << "PC max: " << pc.max() << endl;
	cout << "PC min: " << pc.min() << endl;
	im_file.write(pc, "phasecorrelation.fts");*/
}
