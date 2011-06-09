#ifndef __blend_h
#define __blend_h


#include <image/image.h>
#include <util/mathutil.h>


class Blend
{
  public:
    VISImage<float> _target;
    VISImage<float> _init;
    VISImage<float> _image;
    VISImage<float> _image2;
    VISImage<float> _level_set_change;
    VISImage<float> _level_set_tmp;
    VISImage<float> _level_set_change2;
    VISImage<float> _level_set_tmp2;

// things you need to uppdate alpha
    VISImage<float> _grad_mag, _grad_mag_sq, _total_mag, _total_mag_sq; 
//    VISImage<float> _old_difference, _total_change;

    float _alpha;
    float _dt;
    int _max_x, _max_y;
    
    void calculateChange()
    {
	calculateChange(_alpha);
    }

    void calculateChange(float alpha);
    void calculateChange2(float alpha);
    
    float update() {return(update(FLT_MAX));}
    float update(float max_dt);
    void iterate()
    {
	calculateChange();
	update();
    }


    float iterateMopUp(float max_dt, float sigma)
    {
	VISImage<float> image = _image, image2 = _image2;
	calculateChange(sigma);
	float dt_first = update(max_dt);
	VISImage<float> change = _level_set_change, 
	    change2 = _level_set_change2;
	calculateChange(sigma);
	float dt_second = update(max_dt);

	
	_image = image + (0.5f*dt_first)*change + 
	    (0.5f*dt_second)*_level_set_change;
	_image2 = image2 + (0.5f*dt_first)*change2 + 
	    (0.5f*dt_second)*_level_set_change2;

	return(0.5f*(dt_first + dt_second));
    }

//    float deltaAlpha(float target_difference);

    void reset()
    {
	_image = _init;
	_image2 = _target;
    }
};


#endif
