#include <image.h>
#include <imagefile.h>


#define IMAGE_SIZE 30

main(int argc, char** argv)
{

  // declare an image of type float --- default size is zero
  VISImage<int> image(IMAGE_SIZE,IMAGE_SIZE);

  // declare an object to read files and produce im objects
  VISImageFile im_file;

  int i, j;

  for (j = 0; j < image.height(); j++)
    for (i = 0; i < image.width(); i++)
      {
	cin >> image.poke(i, j);
      }

  im_file.write(image, "image_text.fit");
  cout << "done." << endl;
}

