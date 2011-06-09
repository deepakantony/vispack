
#include "blend.h"


void Blend::calculateChange(float sigma)
{
    VISImage<float> difference_image = _image2 - _image;
    _level_set_tmp2 = 
	_level_set_tmp = 1.0f/((sigma*sigma + difference_image.power(2)).sqrt());

    _level_set_change = difference_image*_level_set_tmp;
    _level_set_change2 = -1.0f*_level_set_change;
	    
//    calculateChange2(sigma);
}



void Blend::calculateChange2(float sigma)
{
    VISImage<float> difference_image = _image - _image2;
    _level_set_tmp2 = 
	1.0f/((sigma*sigma + difference_image.power(2)).sqrt());
    difference_image = difference_image*_level_set_tmp2;

    int w = _image.width(), h = _image.height();
    int i, j, k;
//    old_difference = difference_image;

    _level_set_change2 = difference_image;
}


float Blend::update(float max_dt)
{
    int w = _level_set_change.width(), 
	h = _level_set_change.height();

//    VISImage<float> level_set_tmp2(w, h), level_set_tmp(w, h);

    float grad_before_x, grad_after_x,
	grad_before_y, grad_after_y;


    float grad_x, grad_y, grad_before, grad_after, grad;
// float max_change = 0.0f;
    float x, y, len;
    float this_change;
    int i, j;
    VISImage<float> total_change;


    for (j = 0; j < h; j++) 
	{
	    for (i = 0; i < w; i++) 
		{
		    this_change = _level_set_change.itemAt(i, j);
		    
		    if (i < (w - 1)) 
			grad_before_x = (_image.itemAt(i + 1, j) 
					 - _image.itemAt(i, j));
		    else 
			grad_before_x = (_image.itemAt(i, j) 
					 - _image.itemAt(i - 1, j));
			
		    if (i > 0)
			grad_after_x = (_image.itemAt(i, j) 
					- _image.itemAt(i - 1, j));
		    else 
			grad_after_x = (_image.itemAt(i + 1, j) 
					- _image.itemAt(i, j));

		    if (j < (h - 1)) 
			grad_before_y = (_image.itemAt(i, j + 1) 
					 - _image.itemAt(i, j));
		    else 
			grad_before_y = (_image.itemAt(i, j) 
					 - _image.itemAt(i, j - 1));

		    if (j > 0)
			grad_after_y = (_image.itemAt(i, j)
					- _image.itemAt(i, j - 1) );
		    else
			grad_after_y = (_image.itemAt(i, j + 1)
					- _image.itemAt(i, j) );

		    if (this_change > 0.0f)
			{
			    grad_before = ::VISmax(grad_before_x, 0.0f);
			    grad_after = ::VISmax(-grad_after_x, 0.0f);		       
			    grad_x = grad_after*grad_after + grad_before*grad_before;
			}
		    else if (this_change < 0.0f)
			{
			    grad_before = VISmin(grad_before_x, 0.0f);
			    grad_after = VISmin(-grad_after_x, 0.0f);		       
			    grad_x = grad_after*grad_after + grad_before*grad_before;
			}
		    else 
			grad_x = 0.0f;

			
//	    this is for computing the abs value in order to normalize the
//          curvature term


		    if (this_change > 0.0f)
			{
			    grad_before = ::VISmax(grad_before_y, 0.0f);
			    grad_after = ::VISmax(-grad_after_y, 0.0f);		       
			    grad_y = grad_after*grad_after 
				+ grad_before*grad_before;
			}
		    else if (this_change < 0.0f)
			{
			    grad_before = VISmin(grad_before_y, 0.0f);
			    grad_after = VISmin(-grad_after_y, 0.0f);		       
			    grad_y = grad_after*grad_after 
				+ grad_before*grad_before;
			}
		    else 
			grad_y = 0.0f;

//	    grad = ::VISMAX(grad_x, grad_y);

		    grad = sqrt(grad_x + grad_y);
		    
//		    if (fabs(this_change) > max_change)
//			{
//			    max_change = fabs(this_change);
//			    _max_x = i;
//			    _max_y = j;
//			}
//		    _level_set_change.at(i, j) = grad*this_change;
		    _level_set_tmp.poke(i, j) = 
			grad*(_level_set_tmp.peek(i, j));
		}
	}

// compute the change for the other image
    for (j = 0; j < h; j++) 
	{
	    for (i = 0; i < w; i++) 
		{
		    this_change = _level_set_change2.itemAt(i, j);
		    
		    if (i < (w - 1)) 
			grad_before_x = (_image2.itemAt(i + 1, j) 
					 - _image2.itemAt(i, j));
		    else 
			grad_before_x = (_image2.itemAt(i, j) 
					 - _image2.itemAt(i - 1, j));
			
		    if (i > 0)
			grad_after_x = (_image2.itemAt(i, j) 
					- _image2.itemAt(i - 1, j));
		    else 
			grad_after_x = (_image2.itemAt(i + 1, j) 
					- _image2.itemAt(i, j));

		    if (j < (h - 1)) 
			grad_before_y = (_image2.itemAt(i, j + 1) 
					 - _image2.itemAt(i, j));
		    else 
			grad_before_y = (_image2.itemAt(i, j) 
					 - _image2.itemAt(i, j - 1));
		    
		    if (j > 0)
			grad_after_y = (_image2.itemAt(i, j)
					- _image2.itemAt(i, j - 1) );
		    else
			grad_after_y = (_image2.itemAt(i, j + 1)
					- _image2.itemAt(i, j) );

		    if (this_change > 0.0f)
			{
			    grad_before = ::VISmax(grad_before_x, 0.0f);
			    grad_after = ::VISmax(-grad_after_x, 0.0f);		       
			    grad_x = grad_after*grad_after + grad_before*grad_before;
			}
		    else if (this_change < 0.0f)
			{
			    grad_before = VISmin(grad_before_x, 0.0f);
			    grad_after = VISmin(-grad_after_x, 0.0f);		       
			    grad_x = grad_after*grad_after + grad_before*grad_before;
			}
		    else 
			grad_x = 0.0f;

			
//	    this is for computing the abs value in order to normalize the
//          curvature term


		    if (this_change > 0.0f)
			{
			    grad_before = ::VISmax(grad_before_y, 0.0f);
			    grad_after = ::VISmax(-grad_after_y, 0.0f);		       
			    grad_y = grad_after*grad_after 
				+ grad_before*grad_before;
			}
		    else if (this_change < 0.0f)
			{
			    grad_before = VISmin(grad_before_y, 0.0f);
			    grad_after = VISmin(-grad_after_y, 0.0f);		       
			    grad_y = grad_after*grad_after 
				+ grad_before*grad_before;
			}
		    else 
			grad_y = 0.0f;

//	    grad = ::VISMAX(grad_x, grad_y);

		    grad = sqrt(grad_x + grad_y);
		    
//		    if (fabs(this_change) > max_change)
//			{
//			    max_change = fabs(this_change);
//			    _max_x = i;
//			    _max_y = j;
//			}
//		    _level_set_change2.at(i, j) = grad*this_change;
		    _level_set_tmp2.poke(i, j) = 
			grad*(_level_set_tmp2.peek(i, j));
		}
	}

//    _dt = ::VISmin(0.5f/max_change, max_dt);
    _dt = ::VISmin(0.25f, max_dt);

    float this_change2, norm;
    for (j = 0; j < h; j++) 
	{
	    for (i = 0; i < w; i++) 
		{
		    this_change = _level_set_tmp.itemAt(i, j)
			+ _level_set_tmp2.itemAt(i, j);
		    
		    if ((norm = _dt*this_change) > 1.0f)
			{
			    _level_set_tmp.poke(i, j) /= norm;
			    _level_set_tmp2.poke(i, j) /= norm;
			}
		}
	}
    _level_set_change = (_image2 - _image)*_level_set_tmp;
    _level_set_change2 = (_image - _image2)*_level_set_tmp2;    

//    float max_level_set_change = ::::VISmax(_level_set_tmp.::VISmax(), 
//				       _level_set_tmp2.::VISmax());

    VISImage<float> change2;
//    if (_alpha > 0.0f)
//    _dt = ::VISmin(::VISmin((0.5f/max_change), 1.0f), max_dt);

//    printf("dt is %f \n", _dt);

// must enforce monoticity

//    if (max_level_set_change > 0.0f)
//	_dt = ::VISmin(_dt, 1.0f/max_level_set_change);
    
//    else 
//	_dt = ::VISmin((0.5f/max_change), max_dt);
//    if (max_change > 0.0f)
//	{
    _image +=  
	(total_change = _dt*_level_set_change);
    _image2 +=  
	(change2 = _dt*_level_set_change2);
//	}
//    else 
//	WARN("we got zero change on this update\n");

//    printf("dt is %f \n", _dt);

    return(_dt);

//    return(fsqrt(((total_change + change2)
//		  .power(2)).average()));
}
