#include <image.h>
#include <imagefile.h>
#include <imageRGBA.h>

VISImageRGBA overlayColor(const VISImage<float> &im, 
const VISImage<float> &overlay, rgba color)
{
  if (!im.compareSize(overlay))
    return VISImageRGBA();
  int w = im.width(), h = im.height();
  int i, j;
  VISImageRGBA ret(w, h);
  float im_min = im.min(), im_max = im.max(), factor;
  float over_min = overlay.min(), over_max = overlay.max(), over_factor;
  factor = 192.0f/(im_max - im_min);
  over_factor = 1.0f/(over_max - over_min);
  float over_value, im_value;
  rgba im_color;

  for (j = 0; j < h; j++)
    for (i = 0; i < w; i++)
      {
	over_value = over_factor*(overlay.peek(i, j) - over_min);
	im_value = factor*(im.peek(i, j) - im_min);
	color.at(A) = 255*over_value;
	im_color = rgba(im_value);
	im_color.at(A) = 255*(1 - over_value);
	ret.poke(i, j) = color + im_color;
      }
  return(ret);
}



main(int argc, char** argv)
{

  // declare an image of type float --- default size is zero
  VISImageFile im_file;
  VISImage<float> image;
  VISImage<float> image_mask;
  rgba color(255, 255, 0);

  // declare an object to read files and produce im objects

  if (argc < 3) 
    {
      cout << "Usage: fitstoraw fits_file overlaymask" << endl;
      exit(-1);
    }

  image = VISImage<float>(im_file.read(argv[1]));
  if (!image.isValid())
    {
      cout << "bad fits file: " << argv[1] << endl;
      exit(-1);
    }
  image_mask = VISImage<float>(im_file.read(argv[2]));
  if (!image.isValid())
    {
      cout << "bad fits file: " << argv[1] << endl;
      exit(-1);
    }

  VISImageRGBA im_color = overlayColor(image, image_mask, color);
  im_color.print();

  im_file.write_jpeg(im_color, "coloroverlay.tif", 100);


  cout << "done " << endl;
}

