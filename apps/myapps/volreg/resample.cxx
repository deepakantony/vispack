#include <volume.h>
#include <volumefile.h>
#include <imagefile.h>
#include <param.h>
#include <scan.h>

void insertImage(VISVolume<float> &vol, const VISImage<float> im, int slice)
{
  if (!((vol.height() == im.height())&&(vol.width() == im.width())))
    {
      cout << "bad insert image inputs" << endl;
      exit(-1);
    }
  int i, j;
  int h = im.height(), w = im.width();
  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      vol.poke(i, j, slice) = im.peek(i, j);
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
  char filename[80], prefix[80];


  VISVolume<unsigned short> vol_short;
  VISVolume<float> vol_float, vol_out;
  VISImage<float> im_float;
  int w, h, d;
  int i, j, k;

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
      cout << "about to read file" << endl;
      read_raw(filename, vol_short);
      cout << "done reading file" << endl;
      vol_float = VISVolume<float>(vol_short);
    }
  else
    // if "notValid" is false, that means read didn't work
    {
      cout << "bad volume file specified " << filename << endl;
      exit(-1);
    }

  im_float = vol_float.image(d/2);
  im_file.write(im_float, "resample_in.fit");

  cout << "volshort min and max" << vol_short.min() << " " << vol_short.max() << endl;

  cout << "about to start resampling" << endl;
  int new_w, new_h, new_d;

  int num_scales; 
  if (VPF::set(num_scales, pf["NUM_SCALES"][0]) != VPF::VALID)
    num_scales = 1;
  
  for (k = 0; k < num_scales; k++)
    {
      vol_out = VISVolume<float>(w = ((new_w = w/2) + w%2), 
				 h = ((new_h = h/2) + h%2),
				 d = ((new_d = d/2) + d%2));
      // take averages of subsampled slices and then insert them into the output
      // image
      insertImage(vol_out, 
		  subSampleAverage(0.5f*(vol_float.image(0) + vol_float.image(1))),
		  0);
      for (j = 1; j < new_d; j++)
	{
	  insertImage(vol_out, 
		      subSampleAverage(0.25f*(vol_float.image(2*j-1) 
					      + vol_float.image(2*j+1)) + 
				       0.5f*vol_float.image(2*j)), 
		      j);
	  cout << "done slice " << j <<endl;
	}
      if (d != new_d)
	insertImage(vol_out, 
		    subSampleAverage(0.5f*(vol_float.image(2*d-1) + vol_float.image(2*d - 2))), 
		    d-1);

      cout << "vol float min and max" << vol_out.min() << " " << vol_out.max() << endl;

      if (VPF::set(prefix, pf["OUTPUT_PREFIX"][0]) != VPF::VALID)
	{
	  cout << "bad image file specified " << argv[1] << endl;
	  exit(-1);
	}

      sprintf(filename, "%s_x%dy%dz%d.raw", prefix, w, h, d);
      vol_short = VISVolume<unsigned short>(vol_out);
      cout << "volshort min and max" << vol_short.min() << " " << vol_short.max() << endl;
      write_raw(filename, vol_short);
      //  read_raw(filename, vol_short);
      //  vol_out = VISVolume<float>(vol_short);
      im_float = vol_out.image(new_d/2);
      im_file.write(im_float, "resample_out.fit");
      cout << "vol float min and max" << im_float.min() << " " << im_float.max() << endl;

      cout << "new width " << vol_out.width() << " height " << vol_out.height() << " depth " 
	   << vol_out.depth() << endl;
      vol_float = vol_out;
    }
  cout << "done" << endl;
}

