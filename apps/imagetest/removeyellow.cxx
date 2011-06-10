
// This is a piece of special purpose code to remove the color overlay 

#include <image.h>
#include <imagefile.h>
#include <imagergba.h>


 main(int argc, char** argv)
 {

     // these are the neighbors
     int N[4][2];
     N[0][0] = 1;
     N[0][1] = 0;
     N[1][0] = -1;
     N[1][1] = 0;
     N[2][0] = 0;
     N[2][1] = 1;
     N[3][0] = 0;
     N[3][1] = -1;

   rgba color;
   char* image_name;

   // a typeless base class
   VISIm im;

   // declare an image of type float --- default size is zero
   VISImage<float> image, image_out, image_mask;

   // declare an object to read files and produce im objects
   VISImageFile im_file;

   int i, j, k; 
   char filename[80];

    if ((argc > 1)&&(im = im_file.read(argv[1])).isValid())
     {
       // construct an image of type float (this will cast each pixel)
       image = VISImage<float>(im);
       //       cout << "got an image " << endl;
       //       im.print();
       image_name = argv[1];
     }
   else
     {
       cout << "need input" << endl;
       exit(-1);
     }

    float threshold = 200;

    image_mask = (image > 200.0f);
    im_file.write(image_mask, "image_mask.fts");
    //    cout << "mask min " << image_mask.min() << "  and max " << image_mask.max() << endl;
    int h = image.height(), w = image.width();

    boolean done = false;
    int num_neighbors;
    float total;
    while (!done)
      {
	done = true;
	for (j = 0; j < h; j++)    
	for (i = 0; i < w; i++)
	  {
	    if (image_mask.peek(i, j) > 0.0)
	      {
		//		cout << "got blank " << endl;
		num_neighbors = 0;
		total = 0.0;
		for (k = 0; k < 4; k++)
		  {
		    if (image_mask.peek(i + N[k][0], j + N[k][1]) == 0.0f)
		      {
			num_neighbors++;
			//			cout << "got neighbor " << endl;
			total += image.peek(i + N[k][0], j + N[k][1]);
		      }
		    done = false;
		  }
		if (num_neighbors > 0)
		  {
		    image.poke(i, j) = total/num_neighbors;
		    image_mask.poke(i, j) = 0.0;
		  }
	      }
	  }
      }

    im_file.write(image, "image_out.fts");
    cout << "done." << endl;
 }

