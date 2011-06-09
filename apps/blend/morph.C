#include "blend.h"
#include <image/imagefile.h>
#include <string.h>
#include <util/mathutil.h>
#include <unistd.h>
#include <sys/times.h>

#define IMAGE_WIDTH (200)
#define IMAGE_HEIGHT (200)

float _center_x, _center_y, _radius;

float circle(unsigned x, unsigned y)
{
    return(((x - _center_x)*(x - _center_x) + (y - _center_y)*(y - _center_y)) 
	   < _radius*_radius);
}


#ifdef not_for_now
    VISImage<float> test_image(IMAGE_WIDTH, IMAGE_HEIGHT);
    _radius = IMAGE_WIDTH/3.0f;
    _center_x = IMAGE_WIDTH/2.0f;
    _center_y = IMAGE_HEIGHT/2.0f;
    test_image.evaluate(circle);
    im_file.write(test_image.gauss(1.0f), "circle1.fit");

    _radius = IMAGE_WIDTH/12.0f;
    _center_x = IMAGE_WIDTH/2.0f;
    _center_y = IMAGE_HEIGHT/2.0f;
    test_image.evaluate(circle);
    im_file.write(4.0f*test_image.gauss(1.0f), "circle2.fit");
    
    _radius = IMAGE_WIDTH/4.0f;
    _center_x = 2.0f*IMAGE_WIDTH/5.0f;
    _center_y = 2.0f*IMAGE_HEIGHT/5.0f;
    test_image.evaluate(circle);
    im_file.write(test_image.gauss(1.0f), "circle3.fit");

    _radius = IMAGE_WIDTH/4.0f;
    _center_x = 3.0f*IMAGE_WIDTH/5.0f;
    _center_y = 3.0f*IMAGE_HEIGHT/5.0f;
    test_image.evaluate(circle);
    im_file.write(test_image.gauss(1.0f), "circle4.fit");
#endif

float the_min, the_max;

float sigma = 1.0f;

main(int argc, char** argv)
{

    VISImageFile im_file;
    VISImage<float> init, target;
    float alpha;
    int iterations;
    float k;

    float normalization;

    
    
    if (argc < 4)
	{
	    printf("usage: morph init target iterations\n");
	    exit(-1);
	}

    if (!(init = VISImage<float>(im_file.read(argv[1]))).isValid())
	{
	    printf("failed to read image file\n");
	    exit(-1);
	}

    if (!(target = VISImage<float>(im_file.read(argv[2]))).isValid())
	{
	    printf("failed to read image file\n");
	    exit(-1);
	}

  //   float tau;
//     float total_time;

//     float clk_tck = 1.0f/CLK_TCK;
//     struct tms buf;
//     clock_t start_time;
//     clock_t end_time;
//     float run_time;

//     times(&buf);
//     start_time = buf.tms_utime;
//     printf("start time %d\n", start_time);

    sscanf(argv[3], "%d", &iterations);
//    sscanf(argv[5], "%f", &total_time);
//    sscanf(argv[5], "%f", &tau);
    
    Blend morph;

// need this for the "kerry" image
//    morph._init = 255.0f - init;
//    morph._target = 255.0f - target;
    the_min = ::VISmin(init.min(), target.min());
    the_max = ::VISmax(init.max(), target.max());

    morph._init = (init/(the_max 
			 - the_min))*(0.5f*(init.height() + init.width()));
    morph._target = (target/(the_max - the_min))*(0.5f*(target.height() 
							+ target.width()));

// this says do just the level-set part of the morph
    morph._alpha = -1.0f;
    morph.reset();

//#define num_images (15)

    VISImage<float> *morph_images, *morph_images2;
    

    morph_images = new VISImage<float>[iterations];
    morph_images2 = new VISImage<float>[iterations];

    char str_tmp[80], filename[80];
    
//    int num_between = iterations/num_images;

//    iterations = num_images*num_between;

    float total_difference = sqrt(((morph._target - morph._init)
			      .power(2)).sum());
    float difference_change, this_difference;

//    morph._alpha = 100.0f;

// fix this

//    float time_const = 1.0f;

// get rid of log relationship Ross 7-1-98
//    float delta_alpha = final_alpha/(float)(iterations);

//tau = log(0.1f/time_const);

    float this_time, time_step, time_per_iteration; 
    int i;
//    for (int i = 0; i < iterations; i++)
//	{
//	    total_time += 1.0f/((iterations - i) + 1.0f);
//	    total_time += time_const*
//		exp(tau*(float)i/(float)(iterations - 1.0f));
//	}

//    float time_const = 1.0f/total_time;


//    delta_alpha = 1.0f/(total_time);
//    float delta_alpha = 1.0f;


//    total_time = 0.0f;

    float current_time = 0.0f;
    time_per_iteration = 1.0;
    float change_diff_per_time = -total_difference/(iterations*0.1f);
    float difference_goal, goal_time, last_difference;
    float delta_time;
    float current_difference;
    
    for (i = 0; i < (iterations); i++)
	{
//	    goal_time = 1.0f/((iterations - i) + 1.0f);
//	    goal_time = time_const
//		*exp(tau*(float)i/(float)(iterations - 1.0f));

//	    goal_time = 
//		-1.0f*log(1.0 - (float)(i + 1)/(float)(iterations + 1))

// get rid of log relationship Ross 7-1-98

//		-0.1f*log(1 - (i + 1)*(delta_alpha))
//		-time_const*log(1 - (i + 1)*(delta_alpha))
//		(i + 1)*(delta_alpha)
//		- current_time;



	    
	    this_time = 0.0f;
	    current_difference = this_difference = 
		sqrt(((morph._image - morph._image2).power(2)).sum());

	    difference_goal = this_difference*
		(1.0f - (1.0f/((iterations + 1.0f) - (i + 1.0f))));
	    goal_time = ::VISmax(((difference_goal - this_difference)
			     /change_diff_per_time), 0.0f);

//	    printf("this difference %f and difference goal is %f\n", 
//		   this_difference, difference_goal);
	    
//	    while (this_time < goal_time)
	    while (((this_difference - difference_goal)/current_difference)
		   > 0.0001)
		{
		    last_difference = this_difference;
//		    if (((total_time + this_time) == 0.0f))
//			morph.calculateChange(-1.0f);
// this one clips the difference image
		    morph.calculateChange(0.001f);
//		    else
//			morph.calculateChange
//			    (-log((total_time + this_time)*delta_alpha)
//			     *alpha);
//			    ((1.0f - ((total_time + this_time)*delta_alpha))
//			    *alpha);

//		    printf("total change on iteration %d is %f \n", i+1, 
		    delta_time = 
//			morph.iterateMopUp(goal_time - this_time, 0.0001);
			morph.update(goal_time - this_time);
		    this_time += delta_time;
//		    this_difference = 
//			((morph._image - morph._image2).abs()).sum();
		    this_difference = 
			sqrt(((morph._image - morph._image2).power(2))
			     .sum());
		    change_diff_per_time = 
			(change_diff_per_time + 
			 (this_difference - last_difference)/delta_time)/2.0f;
		    goal_time = this_time +
			::VISmax(((difference_goal - this_difference)
			     /change_diff_per_time), 0.0f);
//		    printf("change per time  %f and diff goal %f\n", 
//			   change_diff_per_time, 
//			   difference_goal - this_difference);		    

//		    printf("%f\t%f\n", 
//			   current_time + this_time,
//			   change_diff_per_time/total_difference
//			   );		    

		    cout <<  "total time " << current_time + this_time
			 << " percent_diff " << change_diff_per_time/total_difference
			 << " this diff " << this_difference << endl;
//			   );		    

//		    printf("goal time is %f and delta time %f\n", 
//			   goal_time, delta_time);
		}

	    cout << "got one " << endl;

	    current_time += this_time;
	    current_difference = this_difference;
	    
//	    printf("total time is %f and goal is %f\n", 
//		   current_time, goal_time);

//	    printf("total time is %f and difference is %f\n", current_time, 
//		   this_difference);



	    
//	    alpha += delta_alpha;

//	    printf("delta alpha is %f \n", 
//		   morph.deltaAlpha((1.0f - 0.01f)*total_difference));
//
//	    exit(-1);
	    
//	    if (((i + 1)%num_between) == 0)
//		{
//		    strcpy(filename, "morph.");
//		    sprintf(str_tmp, "%03d", i+1);
//		    strcat(filename, str_tmp);
//		    strcat(filename, ".fit");
//		    im_file.write(morph._image, 
//				  filename);

//		    morph_images[(i+1)/num_between - 1] = morph._image;
//		    morph_images2[(i+1)/num_between - 1] = morph._image2;

		    morph_images[i] = morph._image;
		    morph_images2[i] = morph._image2;

//		}
	}

//    VISImage<float> morph_difference 
//	= morph._target - morph_images[num_images - 1];


    VISImage<float> difference_image;

    difference_image = 0.5f*(morph_images[iterations - 1] 
			     - morph_images2[iterations - 1]);


    float time;

//     times(&buf);
//     end_time = buf.tms_utime;
//     printf("end time %d\n", end_time);
    
//     printf("Total processing time %3.2f", (float)clk_tck*((float)end_time 
// 					       - (float)start_time));
    

#ifdef not_for_now
    for (i = 0; i < iterations; i++)
	{
//	    time = 
//		-time_const*log(1 - (i + 1)*(delta_alpha));
//		(i + 1)*(delta_alpha);
	    morph_images[i] -= 
		(0.5f*(1.0f - (float)cos(M_PI*i/(iterations - 1))))
		*difference_image;
	    morph_images2[i] += 
		(0.5f*(1.0f - (float)cos(M_PI*i/(iterations - 1))))
		*difference_image;
	}
#endif


    float the_min = FLT_MAX, the_max = -1.0f*FLT_MAX;
    for (i = 0; i < iterations; i++)
	{
	    the_min = ::VISmin(the_min, morph_images[i].min());
	    the_min = ::VISmin(the_min, morph_images2[i].min());
	    the_max = ::VISmax(the_max, morph_images[i].max());
	    the_max = ::VISmax(the_max, morph_images2[i].max());
	}

    the_min = ::VISmin(the_min, (morph._init).min());
    the_min = ::VISmin(the_min, (morph._target).min());
    the_max = ::VISmax(the_max, (morph._init).max());
    the_max = ::VISmax(the_max, (morph._target).max());

    strcpy(filename, "morph.");
    sprintf(str_tmp, "%03d", 0);
    strcat(filename, str_tmp);
    strcat(filename, ".tif");
    im_file.write_tiff(VISImage<byte>(
	255.0f*(morph._init - the_min)/(the_max - the_min)), 
		       filename);

    strcpy(filename, "morph_interp.");
    sprintf(str_tmp, "%03d", 0);
    strcat(filename, str_tmp);
    strcat(filename, ".tif");
    im_file.write_tiff(VISImage<byte>(
	255.0f*(morph._init - the_min)/(the_max - the_min)), 
		  filename);


    float gamma;
    for (i = 0; i < iterations; i++)
	{
	    strcpy(filename, "morph.");
	    sprintf(str_tmp, "%03d", i+1);
	    strcat(filename, str_tmp);
	    strcat(filename, ".tif");
	    im_file.write_tiff(VISImage<byte>(
		255.0f*(morph_images[i] - the_min)/(the_max - the_min)), 
			       filename);

	    strcpy(filename, "morph_interp.");
	    sprintf(str_tmp, "%03d", i+1);
	    strcat(filename, str_tmp);
	    strcat(filename, ".tif");
	    im_file.write_tiff(VISImage<byte>(
		(255.0f*((morph._target*(float)((i+1)/(2.0f*iterations))
			 + 
			   ((morph._init*(float)(((2.0f*iterations) 
						    - (i + 1))
					      /(2.0f*iterations))))
			  - the_min)/(the_max - the_min)))), 
			       filename);
	}

    for (i = 0; i < iterations; i++)
    {
	strcpy(filename, "morph.");
	sprintf(str_tmp, "%03d", i + iterations);
	strcat(filename, str_tmp);
	strcat(filename, ".tif");
	im_file.write_tiff(VISImage<byte>(
	    255.0f*(morph_images2[iterations - (i + 1)] 
		    - the_min)/(the_max - the_min)), 
			   filename);


	strcpy(filename, "morph_interp.");
	sprintf(str_tmp, "%03d", i+iterations);
	strcat(filename, str_tmp);
	strcat(filename, ".tif");
	im_file.write_tiff(VISImage<byte>(
	    (255.0f*((((morph._target*((float)((i+iterations)
					    /(2.0f*iterations))))
		      + 
		    (morph._init*((float)((2.0f*iterations) 
					  - (i + iterations))
				  /(2.0f*iterations)))
		    )
		      - the_min)/(the_max - the_min)))),
			   filename);
	
    }

    strcpy(filename, "morph.");
    sprintf(str_tmp, "%03d", 2*iterations);
    strcat(filename, str_tmp);
    strcat(filename, ".tif");
    im_file.write_tiff(VISImage<byte>(
	255.0f*(morph._target - the_min)/(the_max - the_min)), 
		       filename);

    strcpy(filename, "morph_interp.");
    sprintf(str_tmp, "%03d", 2*iterations);
    strcat(filename, str_tmp);
    strcat(filename, ".tif");
    im_file.write_tiff(VISImage<byte>(
	255.0f*(morph._target - the_min)/(the_max - the_min)), 
		       filename);

    int w = (morph._init).width(), 
	h = (morph._init).height(), j;
    
#ifdef not_for_now
// prints out cross sections of the images
    for (i = 0; i < w; i++)
	{
	    printf("%f\t", (morph._init).itemAt(i, i));
	    for (j = 0; j < iterations; j++)
		printf("%f\t", (morph_images[j]).itemAt(i, i));	  
	    for (j = iterations - 1; j >= 0; j--)
		printf("%f\t", (morph_images2[j]).itemAt(i, i));	  
	    printf("%f\n", (morph._target).itemAt(i, i));
	}
    printf("\n");

    for (i = 0; i < w; i++)
	{
	    printf("%f\t", (morph._init).itemAt(i, h/2));
	    for (j = 0; j < iterations; j++)
		printf("%f\t", (morph_images[j]).itemAt(i, h/2));	  
	    for (j = iterations - 1; j >= 0; j--)
		printf("%f\t", (morph_images2[j]).itemAt(i, h/2));	  
	    printf("%f\n", (morph._target).itemAt(i, i));
	}
    printf("\n");
#endif
//	{
//	    morph_images[i] += ((float)((i+1)*num_between)/(float)iterations)
//		*morph_difference;
//
//	    strcpy(filename, "morph.");
//	    sprintf(str_tmp, "%d", i+1);
//	    strcat(filename, str_tmp);
//	    strcat(filename, ".fit");
//	    im_file.write(morph_images[i], 
//			  filename);
//	}
//    im_file.write(morph._image, "morph.fit");

//     times(&buf);
//     end_time = buf.tms_utime;
//     printf("Total time %3.2f", (float)clk_tck*((float)end_time 
// 					       - (float)start_time)
// 	   /(float)i);

    
}
