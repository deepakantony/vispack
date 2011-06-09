#include <array.h>
#include <param.h>
#include <imagefile.h>
#include <matrix.h>
#include <image.h>
#include "transform.h"

// This is a program that takes as input a single parameter file
// It produces two rotated output images in accordance with the
// parameters given  
// The parameter file has the following paramters
// THETA <float> - rotation angle in radians
// WIDTH <int> - output image width
// HEIGHT <int> - output image height
// SCALE <float> - scale the output image coords (e.g. 2.0 shrinks the
// output) 
// OFFSET <float> <float> - offset the output image (allows you to
// center the warped image in the output)
// IMAGE_FILE <string> - name of image file to be used as input
main(int argc, char** argv)
{
  VISImageFile im_file;
  if (argc < 2)
    {
      cout << "usage: readparams parameterfile" << endl;
      exit(-1);
    }

  float theta, scale;
  VISVector offset(2);
  int w, h;

  // Create a parameter file object with the first argument
  VPF::ParameterFile pf(argv[1]);

  // Get the parameters from the file
  // If one is not given or bad---quit
  if (VPF::set(theta, pf["THETA"][0]) != VPF::VALID)
    {
      cout << "bad rotation theta specified " << argv[1] << endl;
      exit(-1);
    }

  if (VPF::set(w, pf["WIDTH"][0]) != VPF::VALID)
    {
      cout << "bad output width specified " << argv[1] << endl;
      exit(-1);
    }

  if (VPF::set(h, pf["HEIGHT"][0]) != VPF::VALID)
    {
      cout << "bad output width specified " << argv[1] << endl;
      exit(-1);
    }

  if (VPF::set(scale, pf["SCALE"][0]) != VPF::VALID)
    {
      cout << "bad output scale specified " << argv[1] << endl;
      exit(-1);
    }

  if ((VPF::set(offset.poke(0), pf["OFFSET"][0]) != VPF::VALID)||
      (VPF::set(offset.poke(1), pf["OFFSET"][0]) != VPF::VALID))
    {
      cout << "bad output offset specified " << argv[1] << endl;
      exit(-1);
    }

  VISImage<float> im_ret, im_in;
  char filename[80];

  // get the input file name and then read the file (make sure it's valid)
  if (VPF::set(filename, pf["IMAGE_FILE"][0]) != VPF::VALID)
    {
      cout << "bad image file specified " << argv[1] << endl;
      exit(-1);
    }

  // if "notValid" is false, that means read didn't work
  if (!(im_in = VISImage<float>(im_file.read(filename))).isValid())
    {
      cout << "bad image file specified " << filename << endl;
      exit(-1);
    }

  // Create a warp of type "Rotation"
  WarpType< Rotation > rot_transform;

  // set up parameter array
  VISArray<float> params;
  params.poke(0) = theta;
  params.poke(1) = im_in.width()/2;
  params.poke(2) = im_in.height()/2;

  // Give the parameters to the warp object
  rot_transform.setParams(params);
  // create rotated output and save it
  im_ret = rot_transform.resample(w, h, scale, offset, im_in);
  im_file.write(im_ret, "rotation_out.fit");

  // create a warp that consists of two rotations as above
  CascadeWarp cascade_rot(rot_transform, rot_transform);
  im_file.write(cascade_rot.resample
		(w, h, scale, offset, im_in), "rotation2_out.fit");
}
