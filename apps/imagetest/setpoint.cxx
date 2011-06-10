#include <image.h>
#include <imagefile.h>
#include <imagergba.h>


 #define IMAGE_SIZE 200

 void floodFillBoundBox(
			VISImage<float> im_input, VISImage<int> &im_labels, 
			float thresh, int label, 
			int x, int y,
			int &x_lo, int &y_lo, int &x_hi, int &y_hi)
 {
     VISImIndexList list;

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

     int at_x, at_y, next_x, next_y;
     int i;

     x_hi = x_lo = x;
     y_hi = y_lo = y;


 //    printf("x lo %d hi %d y lo %d hi %d\n",
 //	   x_lo, x_hi, y_lo, y_hi);

 //    printf("first new item at x %d y %d with at %f\n", 
 //	   x, y, 
 //	   itemAt(x, y));

     if (!(im_input.itemAt(x, y) >= thresh))
       return;

 //  if it passes the threshold put the starting pixel on the list

     list.appendItem(VISImIndex(x, y));
     x_hi = x_lo + 1;
     y_hi = y_lo + 1;

 //  mark that pixel as "in"
     im_labels.at(x, y) = label;
     list.reset();

     while (list.valid())
	 {
	   at_x = list.atCurrent().a();
	   at_y = list.atCurrent().b();
	   //	    printf("flood fill iteration %d x %d y %d z %d\n", n++, at_x, 
	   //		   at_y );

	   // look at your neighbors and if they are: 1) in bounds, 2) more than the
	   // threshod and 3) not yet visited, put them on the list and mark them as
	   // visited.     
	   for (i = 0; i < 4; i++)
	     {
	       next_x = at_x + N[i][0];
	       next_y = at_y + N[i][1];
	       if (im_input.checkBounds((unsigned)next_x, (unsigned)next_y))
		 {
		   if ((im_labels.itemAt(next_x, next_y) != label)
		       &&(im_input.itemAt(next_x, next_y) >= thresh))
		     {
		       im_labels.at(next_x, next_y) = label;
		       list.appendItem(VISImIndex(next_x, next_y));
		       x_lo = ::VISmin(x_lo, next_x);
		       x_hi = ::VISmax(x_hi, (next_x + 1));
		       y_lo = ::VISmin(y_lo, next_y);
		       y_hi = ::VISmax(y_hi, (next_y + 1));
		     }
		 }
	     }
	   //
	   // remove the guy whose neighbors you have just visited
	   // when the list is empty, you are done.
	   //
	   list.removeCurrent();
	 }
 }


 class Feature
 {
  public:
   int _bbox_w, _bbox_h, _bbox_x, _bbox_y;
   VISImage<float> _mask, _image, _features, _primitive_mask;
   float _center_of_mass_x, _center_of_mass_y;
   float _size, _mass, _confidence, _average_grey;
 };

 VISArray<Feature> makeFeatures(const VISImage<float> &features, 
				const VISImage<float> &im_input, float threshold_lo, 
				float threshold_hi)
{
  int i, j;
  VISArray<Feature> feature_array;
  VISImage<int> labels = features.createToSize();
  labels = -1;

  VISImage<float> im_tmp;
  int w = features.width(), h = features.height();
  int feature_num;
  Feature new_feature;
  int bb_x_lo, bb_x_hi, bb_y_lo, bb_y_hi; 
  int label = 1;
  VISImage<float> im_ramp_x(w, h), im_ramp_y(w, h);
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	im_ramp_x.poke(i, j) = i;
	im_ramp_y.poke(i, j) = j;
      }

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      {
	if ((features.peek(i, j) >= threshold_hi)&&(labels.peek(i, j) == -1))
	  {
	    floodFillBoundBox(features, labels, threshold_lo, label, i, j, 
			      bb_x_lo, bb_y_lo, bb_x_hi, bb_y_hi);

	    cout << "got flood fill at " << i << " " << j << " with bbox " 
		 << bb_x_lo << " " << bb_x_hi << " " 
		 << bb_y_lo << " " << bb_y_hi << endl;
	    bb_x_lo = VISmax(bb_x_lo - 1, 0);
	    bb_y_lo = VISmax(bb_y_lo - 1, 0);
	    bb_x_hi = VISmin(bb_x_hi + 1, w);
	    bb_y_hi = VISmin(bb_y_hi + 1, h);
	    // get the data and the bbox info
	    new_feature._bbox_w = (bb_x_hi - bb_x_lo);
	    new_feature._bbox_h = (bb_y_hi - bb_y_lo);
	    new_feature._bbox_x = bb_x_lo;
	    new_feature._bbox_y = bb_y_lo;
	    cout << "about to do labels" << endl;
	    new_feature._mask = VISImage<float>((labels.getROI
						 (new_feature._bbox_x,
						  new_feature._bbox_y, 
						  new_feature._bbox_w,
						  new_feature._bbox_h) == label));
	    cout << "about to do image" << endl;
	    new_feature._image = im_input.getROI
	      (new_feature._bbox_x,
	       new_feature._bbox_y, 
	       new_feature._bbox_w,
	       new_feature._bbox_h);
	    cout << "about to do features" << endl;
	    new_feature._features = features.getROI
	      (new_feature._bbox_x,
	       new_feature._bbox_y, 
	       new_feature._bbox_w,
	       new_feature._bbox_h);

	    cout << "mask sum " << (new_feature._mask).sum() << endl;
	    cout << "feature sum " << (new_feature._features).sum() << endl;
	    im_tmp = new_feature._features*new_feature._mask;
	    new_feature._mass = im_tmp.sum();
	    if (new_feature._mass > 1.0e-10)
	      {
		new_feature._center_of_mass_x = ((im_ramp_x.getROI(
								   0, 0, 
								   new_feature._bbox_w,
								   new_feature._bbox_h
								   ))*(im_tmp)).sum()/new_feature._mass
		  + new_feature._bbox_x;
		new_feature._center_of_mass_y = ((im_ramp_y.getROI(
								   0, 0, 
								   new_feature._bbox_w,
								   new_feature._bbox_h
								   ))*(im_tmp)).sum()/new_feature._mass
		  + new_feature._bbox_y;
	      }
	    else cout << "WARNING: no features detected" << endl;
	    feature_array.appendItem(new_feature);
	    cout << "Feature " << feature_array.n() << " at " 
		 << new_feature._center_of_mass_x << " " 
		 << new_feature._center_of_mass_y << " with mass " << 
	      new_feature._mass << endl;
	    label++;
	  }
      }
  return(feature_array);
}

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
	 color.at(A) = 255.0*over_value;
	 im_color = rgba(im_value);
	 im_color.at(A) = 255.0*(1 - over_value);
	 ret.poke(i, j) = color + im_color;
       }
   return(ret);
 }

 VISImageRGBA overlayColor(const VISImageRGBA &im, 
 const VISImage<float> &overlay, rgba color)
 {
   if (!im.compareSize(overlay))
     return VISImageRGBA();
   int w = im.width(), h = im.height();
   int i, j;
   VISImageRGBA ret(w, h);
   float over_min = overlay.min(), over_max = overlay.max(), over_factor;
   over_factor = 1.0f/(over_max - over_min);
   float over_value, im_value;
   rgba im_color;

   for (j = 0; j < h; j++)
     for (i = 0; i < w; i++)
       {
	 over_value = over_factor*(overlay.peek(i, j) - over_min);
	 color.at(A) = 255.0*over_value;
	 im_color = im.peek(i, j);
	 im_color.at(A) = 255.0*(1.0 - over_value);
	 ret.poke(i, j) = color + im_color;
       }
   return(ret);
 }


void drawCircle(VISImage<float> &image, float x, float y, float radius, float value)
{
  int i, j;
  for (i = (x - (radius + 2)); i < (x + (radius + 2)); i++)
    for (j = (y - (radius + 2)); j < (y + (radius + 2)); j++)
	      {
		if ((pow((i - x), 2) + pow((j - y), 2)) < radius*radius)
		  if (image.checkBounds(i, j))
		    image.poke(i, j) = value;
	      }
}


 VISImage<float> constructCircleFilter(float radius)
 {
   int width = (int)ceil(2.0f*radius) + 3;
   VISImage<float> circle(width, width);
   int i, j;
   float d;
   float center = width/2.0f;
   for (i = 0; i < width; i++)
     for (j = 0; j < width; j++)
     {
       d = sqrt(power((float)i - center, 2)
		+ power((float)j - center, 2));
       circle.poke(i, j) = d - radius;
     }
   return(VISImage<float>(circle.zeroCrossings()));
 }

 VISImage<float> makeDisk(float radius)
 {
   int width = (int)ceil(2.0f*radius) + 4;
   VISImage<float> disk(width, width);
   int i, j;
   float d;
   float center = width/2.0f - 0.5;
   for (i = 0; i < width; i++)
     for (j = 0; j < width; j++)
     {
       d = sqrt(power((float)i - center, 2)
		+ power((float)j - center, 2));
       disk.poke(i, j) = (float)(d <= radius);
     }
   return(disk);
 }

 VISImage<float> makeCenterSuround(float radius1, float radius2)
 {
   float ratio, sum1, sum2;
   if (!(radius2 > radius1))
     {
       cout << "bad radii in makeCenterSuround" << endl;
       exit (-1);
     }
   VISImage<float> disk1 = makeDisk(radius1);
   VISImage<float> disk2 = makeDisk(radius2);
   int w1 = disk1.width(), w2 = disk2.width();
   VISImage<float> disk_tmp(w2, w2);
   disk_tmp = 0.0f;
   disk_tmp.putROI(disk1, (w2-w1)/2, (w2-w1)/2);
   disk1 = disk_tmp;
   sum1 = disk1.sum();
   sum2 = disk2.sum();
   // ratio of outer ring to inner (inverse)
   ratio = 1.0f/((sum2 - sum1)/sum1);

   disk2 *= ratio;
   disk1 *= (1.0f + ratio);
   disk1 = (disk1 - disk2).gauss(1.0f);
   disk1 /= sqrt((disk1.abs()).sum());
   return(disk1);
 }


 main(int argc, char** argv)
 {
   rgba color;
   char* image_name;

   // a typeless base class
   VISIm im;

   // declare an image of type float --- default size is zero
   VISImage<float> image, image_out;

   // declare an object to read files and produce im objects
   VISImageFile im_file;

   float the_min, the_max;
   int h, w;
   int i, j;
   char filename[80];
   int RADIUS = 35;


    // this is for the void detection.
    int num_disks = 6;
    // set of filters
    VISImage<float> centerS[num_disks];
    for (i = 0; i < num_disks; i++)
      {
        centerS[i] = makeCenterSuround(2*i + 3.0, 2.0*(i + 3.0));
        sprintf(filename, "centerS_%d.fit", i);
        im_file.write(centerS[i], filename);
      }

   if ((argc > 1)&&(im = im_file.read(argv[1])).isValid())
     {
       // construct an image of type float (this will cast each pixel)
       image = VISImage<float>(im);
       cout << "got an image " << endl;
       im.print();
       image_name = argv[1];
     }
   else
     {
       cout << "need input" << endl;
       exit(-1);
     }

   VISImage<float> blurred, voids;
   VISImage<float> edges, circles, circle_filter;
   VISImageRGBA im_color;

    blurred = image.gaussDiffuse(1.5);
    //   blurred = image.gaussDiffuse(0.5);

// get circles from input for debugging.
   if ((argc > 2)&&(im = im_file.read(argv[2])).isValid())
     {
       // construct an image of type float (this will cast each pixel)
       circles = VISImage<float>(im);
     }
   else
     {

    cout << "image height " << image.height() << " and width " << image.width() << endl;

    edges = VISImage<float>(blurred.cannyEdges(2.0));
    im_file.write(edges, "edges_tmp.fts");
    edges = edges.gaussDiffuse(2.0);

    im_file.write(edges, "edges.fits");
    cout << "done writing edges" << endl;

    //    im_file.write(circle_filter = constructCircleFilter(18), "circle3.fits");
    im_file.write(circle_filter = constructCircleFilter(RADIUS), "circle3.fits");
    cout << "about to do convlve" << endl;
    circles = edges.convolve(circle_filter);
    im_file.write(circles, "circles_tmp.fits");
     }

   float max_response;

   max_response = circles.max();
   VISArray<Feature> feature_list;
   feature_list = makeFeatures(circles, image, 0.55*max_response, 
			       0.65*max_response);

   float x, y, value;
   VISImage<float> circle_picts = image.createToSize();
   for (i = 0; i < feature_list.n(); i++)
     {
       x = (feature_list.peek(i))._center_of_mass_x;
       y = (feature_list.peek(i))._center_of_mass_y;
       value = (feature_list.peek(i))._mass;
       if (value > 150.0)
	 drawCircle(circle_picts, x, y, 0.8*RADIUS, 1.0);
     }

   im_file.write(circle_picts, "circle_picts.fits");
   //   circles = (circles > 0.5f*circles.max());
   // good one   circles = (circles > 0.65f*circles.max());


   if ((argc > 3)&&(im = im_file.read(argv[3])).isValid())
     {
       // construct an image of type float (this will cast each pixel)
       voids = VISImage<float>(im);
     }
   else
     {
   // detect voids
    voids = blurred.convolve(centerS[0]);
    cout << "done voids # 0" << endl;
    for (i = 1; i < num_disks; i++)
      {
        voids = voids.max(blurred.convolve(centerS[i]));
        cout << "done voids # " << i << endl;
      }

    voids = voids*circle_picts;
    im_file.write(voids, "voids_tmp.fits");
     }

    voids -= (0.2f)*voids.max();
    voids *= (voids > 0.0f);
    
    blurred = voids;

    voids = blurred.convolve(centerS[0]);
    cout << "done voids # 0" << endl;
    for (i = 1; i < num_disks; i++)
      {
        voids = voids.max(blurred.convolve(centerS[i]));
        cout << "done voids # " << i << endl;
      }

    im_file.write(voids, "voids.fits");
    cout << "done voids" << endl;
 
   im_file.write(circles, "circles2.fits");
   color = rgba(255, 50, 255);
   im_color = overlayColor(image, circles, color);


//    im_file.write(circle_filter = constructCircleFilter(13), "circle1.fits");
//    cout << "about to do convlve" << endl;
//    circles = edges.convolve(circle_filter);
//    circles = (circles > 0.65f*circles.max());
//    im_file.write(circles, "circles1.fits");
//    color = rgba(255, 255, 100);
//    im_color = overlayColor(im_color, circles, color);


   im_file.write_jpeg(im_color, "coloroverlay.tif", 100);

   cout << "done." << endl;
 }

