
#include <imageRGBA.h>
#include <imagefile.h>
#include <param.h>


main(int argc, char** argv)
{
  VISImageFile im_file;
  if (argc < 2)
    {
      cout << "usage: readparams parameterfile" << endl;
      exit(-1);
    }

  // Create a parameter file object with the first argument
  VPF::ParameterFile pf(argv[1]);
  char filename[80];

  VISImage<float> image_1, image_2;
  VISImageRGBA image_color;
  VISImageRGBA colormap(256, 256);
  int i, j;

  // get the input file name and then read the file (make sure it's valid)
  if (argc < 3 )
    {
    if (VPF::set(filename, pf["IMAGE_FILE_1"][0]) != VPF::VALID)
      {
	cout << "bad image file specified " << endl;
	exit(-1);
    }
    }
    else
      sprintf(filename, "%s", argv[2]);
  image_1 = VISImage<float>(im_file.read(filename));
      
  

  // get the input file name and then read the file (make sure it's valid)
  if (argc < 4 )
    {
    if (VPF::set(filename, pf["IMAGE_FILE_2"][0]) != VPF::VALID)
      {
	cout << "bad image file specified " << endl;
	exit(-1);
      }
    }
    else
      sprintf(filename, "%s", argv[3]);


  image_2 = VISImage<float>(im_file.read(filename));

  float the_min, the_max, the_range;

  the_min = image_1.min();
  the_max = image_1.max();
  the_range = the_max - the_min;
  //  the_max -= 0.4*the_range;
  image_1 = (image_1 - the_min)*(1.0f/(the_max - the_min));

//    the_min = VISmin(image_1.min(), image_2.min());
//    the_max = VISmax(image_1.max(), image_2.max());
//    the_range = the_max - the_min;
//    the_min += 0.1*the_range;
//    the_max -= 0.2*the_range;

  the_min = image_2.min();
  the_max = image_2.max();
  image_2 = (image_2 - the_min)*(1.0f/(the_max - the_min));

  float gamma;
  if (VPF::set(gamma, pf["GAMMA_2"][0]) == VPF::VALID)
  // do some gamma correction on 2
    {
      cout << "about to do gamma" << endl;
      image_2 = (gamma*image_2.ln()).exp();
    }

  float angle, radius;
  #define FACTOR (6)
  // Build BY bivariate colormap
  for (i = 0; i < 256; i++)
    for (j = 0; j < 256; j++)
      {
	if ((i == 0)&&(j == 0))
	  colormap.poke(i, j) = rgba(0);
	angle = atan2(i, j);
	radius = exp((1.0f/FACTOR)*log(power((float)i/1.1, FACTOR) + 
				       power((float)j, FACTOR)));
	//	colormap.poke(i, j) = rgba(rint(radius*sin(angle)), 
	//			      rint(radius*sin(angle)), 
	//			      rint(radius*cos(angle)));
	angle = (angle - M_PI/4.0f)/(M_PI/4.0f);
	if (angle != 0.0f)
	  angle *= angle*angle/fabs(angle);
	angle = angle*M_PI/4.0f + M_PI/4.0f;
	colormap.poke(j, i) = rgba(rint(radius*sin(angle)), 
				   //				   0,
				   rint(radius*sin(angle)), 
				   rint(radius*cos(angle)));
	
      }
  int w = image_1.width(), h = image_1.height();
  image_color = VISImageRGBA(w, h);
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      image_color.poke(i, j) 
	= colormap.peek(VISmax(VISmin((int)rint(255.0f*(image_2.peek(i, j))), 255), 0), 
			VISmax(VISmin((int)rint(255.0f*(image_1.peek(i, j))), 255), 0));
  im_file.write_jpeg(image_color, "WB-BB_bivariate.jpg", 100);
  im_file.write_jpeg(colormap, "colormap.jpg", 100);

  // make checkboard image
  int pixels_per_check;
  VISImageRGBA im_check(w, h);
  if (VPF::set(pixels_per_check, pf["PIXELS_PER_CHECK"][0]) == VPF::VALID)
    {
      for (i = 0; i < w; i++)
	for (j = 0; j < h; j++)
	  if (((i/pixels_per_check)%2 == 0)==((j/pixels_per_check)%2 == 0))
	    im_check.poke(i, j) 
	      = rgba(0,
		     (byte)rint(255.0f*VISmax(VISmin(image_1.peek(i, j), 1.0f), 0.0f)), 
		     (byte)rint(255.0f*VISmax(VISmin(image_1.peek(i, j), 1.0f), 0.0f)));
	  else
	    im_check.poke(i, j) 
	      = rgba((byte)rint(255.0f*VISmax(VISmin(image_2.peek(i, j), 1.0f), 0.0f)),
		     0, 0);
      if (argc > 4)
	sprintf(filename, "image_check_%s.jpg", argv[4]);
      else
	sprintf(filename, "image_check.jpg");
      im_file.write_jpeg(im_check, filename, 100);

      pixels_per_check /= 4;
      for (i = 0; i < w; i++)
	for (j = 0; j < h; j++)
	  if (((i/pixels_per_check)%2 == 0)==((j/pixels_per_check)%2 == 0))
	    im_check.poke(i, j) 
	      = rgba(0,
		     (byte)rint(255.0f*VISmax(VISmin(image_1.peek(i, j), 1.0f), 0.0f)), 
		     (byte)rint(255.0f*VISmax(VISmin(image_1.peek(i, j), 1.0f), 0.0f)));
	  else
	    im_check.poke(i, j) 
	      = rgba((byte)rint(255.0f*VISmax(VISmin(image_2.peek(i, j), 1.0f), 0.0f)),
		     0, 0);
      im_file.write_jpeg(im_check, "image_check2.jpg", 100);
    }

//    i = 635; j = 222;
//    cout << "color at " << i << " " 
//         << j << " is " << image_1.peek(i, j) 
//         << " " << image_2.peek(i, j)  << endl;

  VISImage<float> histogram(128, 128);
  // build the joint histogram
  // Build BY bivariate colormap
  for (i = 0; i < 128; i++)
    for (j = 0; j < 128; j++)
      {
	histogram.poke(i, j) = 0.0f;
      }

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      histogram.poke(VISmax(VISmin((int)rint(127.0f*(image_1.peek(i, j))), 127), 0), 
		     VISmax(VISmin((int)rint(127.0f*(image_2.peek(i, j))), 127), 0))
	+= 1;

  histogram =   (1.0f + histogram.gaussDiffuse(1.0)).ln();
  
  im_file.write(histogram, "histogram.fit");

//    for (i = 0; i < 128; i++)
//        {
//  	for (j = 0; j < 128; j++)
//  	  cout << log(1.0 + histogram.peek(i, j)) << " ";
//  	cout << endl;
//        }

}


