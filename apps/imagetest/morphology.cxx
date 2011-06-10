#include <image.h>
#include <imagefile.h>


VISImage<byte> createCircleMask(int radius)
{
  if (radius < 0)
    return(VISImage<byte>());
  int w, h;
  w = h = radius*2 + 1;
  int i, j;
  int center_x = w/2, center_y = h/2;
  VISImage<byte> ret(w, h);
  for (j = 0; j < h; j++)
    for (i = 0; i < w; i++)
      {
	if ((power(i - center_x, 2) + power(j - center_y, 2) <= radius*radius))
	  ret.poke(i, j) = 1;
	else 
	  ret.poke(i, j) = 0;
      }
  return(ret);
}

VISImage<byte> dilate(const VISImage<byte> &mask, const VISImage<byte> &image)
{
  int w = image.width(), h = image.height();
  VISImage<byte> ret(w, h);
  ret = (byte)0;
  int mask_w = mask.width()/2, mask_h = mask.height()/2;
  
  int i, j, ii, jj, new_i, new_j;
  boolean condition;

  for (j = 0; j < h; j++)  
    for (i = 0; i < w; i++)
      {
	condition = false;
	for (jj = -mask_h; jj <= mask_h; jj++)  
	  for (ii = -mask_w; ii <= mask_w; ii++)  
	    {
	      new_i = i + ii; new_j = j + jj;
	      if (image.checkBounds(new_i, new_j))
		{
		  if (mask.peek(ii + mask_w, jj + mask_h)&&image.peek(new_i, new_j))
		    condition = true;
		}
	    }
	ret.poke(i, j) = (byte)condition;
      }
  return(ret);
}

VISImage<byte> erode(const VISImage<byte> &mask, const VISImage<byte> &image)
{
  int w = image.width(), h = image.height();
  VISImage<byte> ret(w, h);
  ret = (byte)0;
  int mask_w = mask.width()/2, mask_h = mask.height()/2;
  
  int i, j, ii, jj, new_i, new_j;
  boolean condition;

  for (j = 0; j < h; j++)  
    for (i = 0; i < w; i++)
      {
	condition = true;
	for (jj = -mask_h; jj <= mask_h; jj++)  
	  for (ii = -mask_w; ii <= mask_w; ii++)  
	    {
	      new_i = i + ii; new_j = j + jj;
	      if (image.checkBounds(new_i, new_j))
		{
		  if (mask.peek(ii + mask_w, jj + mask_h)&&(!image.peek(new_i, new_j)))
		    condition = false;
		}
	    }
	ret.poke(i, j) = (byte)condition;
      }
  return(ret);
}

VISImage<byte> close(const VISImage<byte> &mask, const VISImage<byte> &image)
{
  return(erode(mask, dilate(mask, image)));
}

VISImage<byte> open(const VISImage<byte> &mask, const VISImage<byte> &image)
{
  return(dilate(mask, erode(mask, image)));
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

  VISImage<byte> segmentation;
  float threshold = 165.0;

  VISImage<byte> mask, im_tmp;

    #ifdef __isotropic_filter
  segmentation = image > threshold;
  im_file.write(segmentation, "outfile.fts");
  mask = createCircleMask(1);
  mask.printData();
  segmentation = open(mask, segmentation);
  im_file.write(segmentation, "outfile1.fts");
  segmentation = close(mask, segmentation);
  im_file.write(segmentation, "outfile2.fts");
  mask = createCircleMask(3);
  mask.printData();
  segmentation = open(mask, segmentation);
  im_file.write(segmentation, "outfile3.fts");
  segmentation = close(mask, segmentation);
  im_file.write(segmentation, "outfile4.fts");
    #endif
  

    //#ifdef orientedfilter

  threshold = 121.0;
  segmentation = image < threshold;
  im_file.write(segmentation, "outfile.fts");
  // use a horizontal structuring element...
  mask = VISImage<byte>(3, 1);
  mask = (byte)1;
  im_tmp = open(mask, segmentation);
  im_file.write(im_tmp, "outfile1.fts");

  mask = VISImage<byte>(1, 3);
  mask = (byte)1;
  segmentation = open(mask, segmentation);
  im_file.write(segmentation, "outfile2.fts");

  segmentation += im_tmp;
  im_file.write(segmentation, "outfile3.fts");

  //  mask = createCircleMask(1);
  //  segmentation = open(mask, segmentation);
  //  im_file.write(segmentation, "outfile4.fts");

  mask = VISImage<byte>(3, 1);
  mask = (byte)1;
  im_tmp = close(mask, segmentation);
  im_file.write(im_tmp, "outfile5.fts");

  mask = VISImage<byte>(1, 3);
  mask = (byte)1;
  im_tmp = close(mask, im_tmp);
  im_file.write(im_tmp, "outfile6.fts");

  //  #endif
  
}

