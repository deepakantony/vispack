/* 
   viewer.cxx
   author: Sarah Geneser
   written: 06-27-03
   last updated: 06-27-03

   inputs a parameter file that specifies two images and a transform and creates 
   a four part image for the purpose of visualizing how good a registration the 
   transform provides.

   upper left: reference image
   lower left: check of reference and transformed test   
   upper right: transformed test image   
   lower right: difference of reference and transformed test

*/



#include <math.h>
#include "image/image.h"
#include "image/imagefile.h"
#include "util/mathutil.h"
#include <string>
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <volume.h>
#include <volumefile.h>
#include <imagefile.h>
#include <param.h>
#include "transform.h"

#define DEBUG
#define SLICE 0


string getFileName(string file)
{
  int begin = file.find_last_of("/");
  string fileName = file.substr(begin, file.size());
  return fileName;
}


// A helper function which reads, from a ParameterFile object, a
// vector of the form (name elements ....data... )
// Returns an invalid vector if it fails (prints something)
VISVector readVector(VPF::ParameterFile pf, const char* name, boolean &good)
{
  int n, j;
  float float_tmp;
  int int_tmp;
  // declare the return object
  VISVector vector;
  // get the number of elements
  good=true;
  if (VPF::set(n, pf[name][0]) != VPF::VALID)
    good=false;

  if (!good)
  {
    cout << "ImplicitModel::bad vector " << name << endl;
    return VISVector();  
  }

  // create the return object
  vector = VISVector(n);
  // read the data
  for (j = 0; j < n; j++)
    if (VPF::set(float_tmp, pf[name][j+2]) == VPF::VALID)
      vector.poke(j) = float_tmp;
    else 
      if (VPF::set(int_tmp, pf[name][j+2]) == VPF::VALID)
	vector.poke(j) = (float)int_tmp;
      else
	good=false;
  
  if (!good)
    {
      cout << "ImplicitModel::bad vector " << name << endl;
      return VISVector();
    }
  return(vector);
}

//function that adds a grid of points to an image in a specified
//region of interest
VISImage<float> addGrid(VISImage<float> im, VISVector lower, VISVector upper,\
			int dots, float color)
{
  int checkWidth = (int)floor((upper.peek(0)-lower.peek(0))/dots),\
    checkHeight = (int)floor((upper.peek(1)-lower.peek(1))/dots);
  int xStart = (int)lower.peek(0),  xEnd = (int) upper.peek(0),\
    yStart = (int)lower.peek(1),  yEnd = (int)upper.peek(1);
  
  for (int i=xStart; i<=xEnd; i=checkWidth+i)
  {
    for (int j=yStart; j<=yEnd; j=checkHeight+j)
    {
      im.poke(i,j) = color;
    }
  }
  return im;
}

//function that inputs a reference image and a test image 
//that has already been transformed and creates a checkered
//combination of the two images with <checks> number of checkers
VISImage<float> makeChecks(VISImage<float> ref, VISImage<float> test, int checks)
{
  //check that the images are the same size
  if ( !( (test.width() == ref.width()) && (test.height() == ref.height()) ) )
    cout << "ERROR: images are assumed to be the same size" << endl;
    
  int w = ref.width(), h = ref.height();
  VISImage<float> ret(w, h);
  ret = 1.0f;
  
  // determine the size of each check
  int checkWidth = (int)floor((float)w/checks), checkHeight = (int)floor((float)h/checks);
  VISImage<float> checker(checkWidth, checkHeight,1);
  
  bool which = true;
  for (int i=0; i<checks; i++)
  {
    for (int j=0; j<checks; j++)
    {
      w = i*checkWidth; h = j*checkHeight;
      //grab a check from the test or reference image
      if (which)
	checker = test.getROI(w, h, checkWidth, checkHeight);
      else 
	checker = ref.getROI(w, h, checkWidth, checkHeight);
      //put the check into the output image
      ret.putROI(checker, w, h);
      //switch between test and reference
      which = (!which);
    }
  }
  return ret;
}

//map the pixel values between 0 and 1
VISImage<float> normalize(VISImage<float> im)
{
  im -= im.min();
  im /= im.max();
  return im;
}
    
//function that inputs a reference image and a test image 
//that has already been transformed and creates the difference
//image between the two
VISImage<float> makeDiff(VISImage<float> ref, VISImage<float> test)
{
  //check that the images are the same size
  if ( !( (test.width() == ref.width()) && (test.height() == ref.height()) ) )
    cout << "ERROR: images are assumed to be the same size" << endl;
  
  int w = ref.width(), h = ref.height();
  VISImage<float> ret(w, h);
  
  ret = test -ref;
  return ret;
}

main (int argc, char **argv)
{

  if (argc < 2) 
    cout << "ERROR: must specify a parameter file" << endl;

  VISVolume<float> testVol, refVol;
  VISImage<float> test, ref, check, difference;
  VISVolumeFile vol_file;
  VISImageFile im_file;
  VISIm im;
  int slice = 0;
  
  //declare VPF
  
  VPF::ParameterFile pf(argv[1]);

#ifdef DEBUG
  cout << "Reading in the parameter file..." << endl;
#endif
  
  //read in the reference image
  char refName[100];
  if (VPF::set(refName, pf["REF_FILE"][0]) == VPF::VALID)
  {
    refVol = vol_file.read_float(refName);
    ref = refVol.image(slice);
  }
  else
    cout << "Warning! invalid reference volume filename specified" << endl;
#ifdef DEBUG
  cout << "Read in the reference volume" << endl;
#endif
  
  //read in the test image
  char testName[500];
  if (VPF::set(testName, pf["TEST_FILE"][0]) == VPF::VALID)
  {
    testVol = vol_file.read_float(testName);
    test = testVol.image(slice);
  }
   else
     cout << "Warning! invalid test volume filename specified" << endl;
#ifdef DEBUG
  cout << "Read in the test volume" << endl;
#endif

  //read in the output directory
  char saveDir[500];
  if (!(VPF::set(saveDir, pf["OUTPUT_DIR"][0]) == VPF::VALID))
    cout << "Warning! improperly specified output directory" << endl;
  
  //read in the initial transform
  //EulerRigidTransform *transform = new EulerRigidTransform();
  AffineTransform *transform = new AffineTransform();
  //TranslateTransform *transform = new TranslateTransform();
  
  boolean good = true;
  VISVector params = readVector(pf, "TRANS", good);
  if (good)
    transform->setParams(params);
  else
    cout << "Warning! invalid transform specified" << endl;
#ifdef DEBUG
  cout << "Read in the transform" << endl;
#endif
 
  int checks;
  //read in the number of checks
  if (!(VPF::set(checks, pf["CHECKS"][0]) == VPF::VALID))
    cout << "Warning! invalid number of checks specified" << endl;
#ifdef DEBUG
  cout << "Number of checks specified is: " << checks << endl;
#endif

  int dots;
  //read in the number of dots
  if (!(VPF::set(dots, pf["DOTS"][0]) == VPF::VALID))
    cout << "Warning! invalid number of dots specified" << endl;
#ifdef DEBUG
  cout << "Number of dots specified is: " << dots << endl;
#endif

  boolean alsoGood;
  //Read in the lower point for the region of interest
  VISVector lower = readVector(pf, "LOWER_ROI", good);
  //Read in the upper point for the region of interest
  VISVector upper = readVector(pf, "UPPER_ROI", alsoGood);
  if(!alsoGood || ! good)
  {
    cout << "Warning! ROI incorrectly specified, defaulting"
	 << " to the entire image" << endl;
    lower = VISVector(2); lower = 0.0f;
    upper = VISVector(2); upper.poke(0) = ref.width(); upper.poke(1) = ref.height();
  }
#ifdef DEBUG
  cout << "Read in the region of interest" << endl;
#endif

  //check that test and ref are the same size
  if ( !( (test.width() == ref.width()) && (test.height() == ref.height()) ) )
    cout << "ERROR: images are assumed to be the same size" << endl;

#ifdef DEBUG
  cout << "Finished reading in the parameter file." << endl;
#endif
  
  int width = test.width(), height = test.height();
  VISImage<float> output(2*width, 2*height);
  
  //transform the test image via the specified transform
  testVol = transform->transformVolume(testVol);
  test = testVol.image(slice);
  test = normalize(test);
  
  //write out the test transform
  char transString[500];
  char viewerString[500];
  string fileName=getFileName(testName);
  sprintf(viewerString, "%s/viewer.fits", saveDir, fileName.c_str());
  sprintf(transString, "%s%s_transformed.fits", saveDir, fileName.c_str());
  im_file.write_fits(test, transString);
  cout << "Finished writing out transformed image to: " << transString << endl;

  ref = normalize(ref);
  
  //create the checkered image
  ////////////////////////////
  check = makeChecks(ref, test, checks);
  //  check = normalize(check);
  //and add that image to the output
  output.putROI(check, 0, height);
#ifdef DEBUG
  cout << "Finished with checkered image" << endl;
#endif

  //create the difference image
  /////////////////////////////
  difference = makeDiff(ref, test);
  difference = normalize(difference);
  //add that image to the output
  output.putROI(difference, width, height);
#ifdef DEBUG
  cout << "Finished with difference image" << endl;
#endif

  float white = ref.max();

  //add the grid in the ROI to the reference image
  ////////////////////////////////////////////////
  ref = addGrid(ref, lower, upper, dots, white);
  ref = normalize(ref);
  //put the reference image into the output
  output.putROI(ref,0,0);
#ifdef DEBUG
  cout << "Finished adding grid to reference image" << endl;
#endif

  white = test.max();
  //add the grid in the ROI to the test image
  ///////////////////////////////////////////
  test = addGrid(test, lower, upper, dots, white);
  test = normalize(test);
  //and add that image to the output
  output.putROI(test, width, 0);
#ifdef DEBUG
  cout << "Finished adding grid to test image" << endl;
#endif
  
  //write out the registration viewer
  im_file.write_fits(output, viewerString);
  cout << "finished writing viewer to: " << viewerString << endl;

  cout << "done." << endl;
}

