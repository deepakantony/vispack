#include <image.h>
#include <imagefile.h>

main(int argc, char** argv)
{
  VISIm im;
  VISImage<float> image1, image2, image3;
  VISImageFile im_file;

  if ((argc > 1)&&(im = im_file.read(argv[1])).isValid())
    {
      // construct an image of type float (this will cast each pixel)
      image1 = VISImage<float>(im);
    }
  else
    {
      cout << "Bad input file 1" << endl; 
      exit(-1);
    }
  if ((argc > 2)&&(im = im_file.read(argv[2])).isValid())
    {
      // construct an image of type float (this will cast each pixel)
      image2 = VISImage<float>(im);
    }
  else
    {
      cout << "Bad input file 2" << endl; 
      exit(-1);
    }

  //  image1 = image1.gaussDiffuse(3.0);
  //  image2 = image2.gaussDiffuse(3.0);
  //  image3 = (image1 - image2).abs();
  im_file.write(image1 + noiseGauss(image1.width(), image1.height(), 50.0f), "image_diff.fts");
  cout << "done" << endl;
}
