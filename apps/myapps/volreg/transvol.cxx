
#include "transform.h"
#include <volumefile.h>
#include <param.h>


void volConvertBL(VISVolume<unsigned short> &vol)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  char *buf_in = (char*)(vol.repRef())->bufferRef();
  char char_tmp;
  cout << "about to do convert " <<  w*h*d << endl;
  cout << "about to do convert " <<  w << " " << h << " " << d << endl;
  for (int k = 0; k < w*h*d; k++)
    {
      if (k%(w*h) == 0)
	cout << "about to do convert " << k << endl;
      char_tmp = buf_in[2*k+1];
      buf_in[2*k+1] = buf_in[2*k];
      buf_in[2*k] = char_tmp;
    }
}

main(int argc, char** argv)
{
  VISVolumeFile vol_file;
  VISImageFile im_file;
  if (argc < 2)
    {
      cout << "usage: readparams parameterfile" << endl;
      exit(-1);
    }

  // Create a parameter file object with the first argument
  VPF::ParameterFile pf(argv[1]);
  char filename[80];


  VISVolume<unsigned short> vol_short;
  VISVolume<int> vol_int;
  VISVolume<float> vol_in, vol_reg, vol_ret;
  int w, h, d;

  // get the input file name and then read the file (make sure it's valid)
  if (VPF::set(filename, pf["VOLUME_FILE"][0]) != VPF::VALID)
    {
      cout << "bad image file specified " << argv[1] << endl;
      exit(-1);
    }
   if ((VPF::set(w, pf["INPUT_SIZE"][0]) == VPF::VALID)
       &&(VPF::set(h, pf["INPUT_SIZE"][1]) == VPF::VALID)
       &&(VPF::set(d, pf["INPUT_SIZE"][2]) == VPF::VALID))
     {
       vol_short = VISVolume<unsigned short>(w, h, d);
       read_raw(filename, vol_short);
       //      vol_short /= (unsigned short)255;
       //      vol_int = VISVolume<int>(vol_short);
       //      cout << vol_int.min() << " " << vol_int.max() << endl;
       //       volConvertBL(vol_short);
       vol_in = VISVolume<float>(vol_short);
       //      vol_in /= (255.0f*255.0f);
     }
   else
  // if "notValid" is false, that means read didn't work
  if (!(vol_in = VISVolume<float>(vol_file.read(filename))).isValid())
    {
      cout << "bad volume file specified " << filename << endl;
      exit(-1);
    }

   //  cout << vol_in.min() << " " << vol_in.max() << endl;

  VISImage<float> im_float = vol_in.image(vol_in.depth()/2);
  //  cout << im_float.min() << " " << im_float.max() << endl;

  //  vol_file.write_float(vol_in, "wb_float.vol");
  im_file.write(im_float, "vol_in_slice.fits");
  

  // get the input file name and then read the file (make sure it's valid)
  if ((VPF::set(w, pf["OUTPUT_SIZE"][0]) != VPF::VALID)
      ||(VPF::set(h, pf["OUTPUT_SIZE"][1]) != VPF::VALID)
      ||(VPF::set(d, pf["OUTPUT_SIZE"][2]) != VPF::VALID))
    {
      cout << "bad volume size specified " << filename << endl;
      exit(-1);
    }

  VISArray<float> affine_params;

  int i;
  float float_tmp;

  for (i = 0; i < 12; i++)
    if (VPF::set(float_tmp, pf["A_TRANSFORM"][i]) != VPF::VALID)
      {
	cout << "bad volume file specified " << argv[1] << endl;
	exit(-1);
      }
    else
      {
	affine_params.poke(i) = float_tmp;
      }

  VISVector origin_source((float)vol_in.width()/2.0f, 
			  (float)vol_in.height()/2.0f, 
			  (float)vol_in.depth()/2.0f),
    origin_target((float)w/2.0f, (float)h/2.0f, (float)d/2.0f);
  // Create a warp of type "Rotation"
  WarpType< Affine > affine_warp;
  affine_warp.setParams(affine_params);
  vol_ret = affine_warp.resample(w, h, d, 1.0f, origin_target, origin_source, vol_in);
  //  vol_file.write_float(vol_ret, "affine_out.vol");
  vol_short = VISVolume<unsigned short>(vol_ret);
  write_raw("registered_wb_x724y512z120_ushort_le.raw", vol_short);
  VISImage<float> im_affine;
  im_file.write((im_affine = vol_ret.image(d/2)), "WB.fits");

  // get the input file name and then read the file (make sure it's valid)
  if (VPF::set(filename, pf["COMPARE_FILE"][0]) != VPF::VALID)
    {
      cout << "bad image file specified " << argv[1] << endl;
      exit(-1);
    }
  if ((VPF::set(w, pf["COMPARE_SIZE"][0]) == VPF::VALID)
      &&(VPF::set(h, pf["COMPARE_SIZE"][1]) == VPF::VALID)
      &&(VPF::set(d, pf["COMPARE_SIZE"][2]) == VPF::VALID))
     {
       vol_short = VISVolume<unsigned short>(w, h, d);
       read_raw(filename, vol_short);
       //      vol_short /= (unsigned short)255;
       //      vol_int = VISVolume<int>(vol_short);
       //      cout << vol_int.min() << " " << vol_int.max() << endl;
       //       volConvertBL(vol_short);
       vol_in = VISVolume<float>(vol_short);
       //      vol_in /= (255.0f*255.0f);
     }
  else
  // if "notValid" is false, that means read didn't work
    if (!(vol_in = VISVolume<float>(vol_file.read(filename))).isValid())
      {
	cout << "bad volume file specified " << filename << endl;
	exit(-1);
      }
  im_file.write((vol_in.image(d/2)), "BB.fits");
  im_file.write((im_affine - vol_in.image(d/2)), "compare_out.fits");
}



