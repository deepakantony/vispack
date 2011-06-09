#include "transform.h"
#include <volumefile.h>
#include <imagefile.h>
#include <param.h>

main(int argc, char** argv)
{
  VISVolume<unsigned short> vol_short;
  VISVolume<float> vol_float;
  VISImageFile im_file;
  VPF::ParameterFile pf(argv[1]);
  
  char filename[80], prefix[80];
  int w, h, d;
  int i;
  if (VPF::set(filename, pf["INPUT_FILE"][0]) != VPF::VALID)
    {
      cout << "bad image file specified " << argv[1] << endl;
      exit(-1);
    }

   if (!((VPF::set(w, pf["INPUT_SIZE"][0]) == VPF::VALID)
       &&(VPF::set(h, pf["INPUT_SIZE"][1]) == VPF::VALID)
	 &&(VPF::set(d, pf["INPUT_SIZE"][2]) == VPF::VALID)))
      cout << "bad image specified " << filename << endl;
  // if "notValid" is false, that means read didn't work
       vol_short = VISVolume<unsigned short>(w, h, d);
       read_raw(filename, vol_short);
       vol_float = VISVolume<float>(vol_short);
       cout << "input min " << vol_float.min() << " and max " << vol_float.max() << endl;

       if (VPF::set(prefix, pf["OUTPUT_PREFIX"][0]) != VPF::VALID)
	 {
	   cout << "bad output file prefix" << endl;
	   exit(-1);
	 }
       
       for (i = 0; i < d; i++)
	 {
	   sprintf(filename, "%s%03d.fit", prefix, i);
	   im_file.write(vol_float.image(i), filename);
	 }
       
}
