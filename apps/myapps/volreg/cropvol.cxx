
#include "transform.h"
#include <volumefile.h>
#include <param.h>

void volConvertBL(VISVolume<unsigned short> &vol)
{
  int w = vol.width(), h = vol.height(), d = vol.depth();
  char *buf_in = (char*)(vol.repRef())->bufferRef();
  char char_tmp;
  //  cout << "about to do convert " <<  w*h*d << endl;
  //  cout << "about to do convert " <<  w << " " << h << " " << d << endl;
  for (int k = 0; k < w*h*d; k++)
    {
      //      if (k%(w*h) == 0)
      //	cout << "about to do convert " << k << endl;
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


  VISVolume<unsigned short> vol_short, vol_short2;
  VISVolume<float> vol_float;
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
       read_raw(filename, vol_short);
       //      vol_short /= (unsigned short)255;
       //      cout << vol_int.min() << " " << vol_int.max() << endl;
       //       volConvertBL(vol_short);
       vol_short = vol_short.getROI(0, 150, 0, vol_short.width() , 
				    vol_short.height() - 300, vol_short.depth());
       if (VPF::set(filename, pf["OUTPUT_FILE"][0]) != VPF::VALID)
	 {
	   cout << "bad image file specified " << argv[1] << endl;
	   exit(-1);
	 }

       // swap x and y because ???????
       vol_short2 = vol_short;
       h = h - 300;
       vol_short = VISVolume<unsigned short>(h, w, d);
       for (i = 0; i < w; i++)
	 for (j = 0; j < h; j++)
	   for (k = 0; k < d; k++)
	     vol_short.poke(j, i, k) = vol_short2.peek(i, j, k);

       vol_float = VISVolume<float>(vol_short);
       im_file.write(vol_float.image(d/2), "cropped_out.fits");
       write_raw(filename, vol_short);
     }
   else
  // if "notValid" is false, that means read didn't work
     {
       cout << "bad volume file specified " << filename << endl;
       exit(-1);
     }

   
}



