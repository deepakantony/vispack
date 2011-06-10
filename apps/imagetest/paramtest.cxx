#include <array.h>
#include <param.h>
#include <imagefile.h>
#include <matrix.h>
#include <image.h>

main(int argc, char** argv)
{
  int w, h;
  float theta;
  cout<<"width, height initialized to "<<w<<"  "<<h<<endl;

  // Create a parameter file object with the first argument
  VPF::ParameterFile pf(argv[1]);

  // Get the parameters from the file
  // If one is not given or bad---quit
  if (VPF::set(theta, pf["THETA"][0]) != VPF::VALID)
    {
      cout << "bad rotation theta specified " << argv[1] << endl;
      cout << "read theta: "<<theta<<endl;
    }
  else cout << "got good theta " << theta << endl;

  if (VPF::set(h, pf["HEIGHT"][0]) != VPF::VALID)
    {
      cout << "bad output height specified " << argv[1] << endl;
      cout << "read height: "<<w<<endl;
    }
  else cout << "got good height " << h << endl;

  if (VPF::set(w, pf["WIDTH"][0]) != VPF::VALID)
    {
      cout << "bad output width specified " << argv[1] << endl;
      cout << "read width: "<<w<<endl;
    }
  else cout << "got good width " << w << endl;

}

