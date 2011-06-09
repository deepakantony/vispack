/*
  MIregister.cxx
  Sarah Geneser
  05-21-03
  uses Mutual Information (Viola and Wells) to register 2D images
  requires a parameter file (param.txt)
  and outputs a parameter file to generate a veiwer of the best 
  registration found.  (viewer viewerParams.txt)
*/


#include "util/mathutil.h"
#include <matrix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <image.h>
#include <volume.h>
#include <volumefile.h>
#include <imagefile.h>

#include "registrator.h"

#include <param.h>


#define DEBUG 1
//#define WARNING 1
#define TWOPI 6.283185307179586
#define eps .000000001


// A helper function which reads, from a ParameterFile object, a
// matrix of the form (name rows cols ....data in row major order...)
// Returns an invalid matrix if it fails (prints something)
VISMatrix readMatrix(VPF::ParameterFile pf, const char* name, boolean &good)
{
  int r, c, j;
  float float_tmp;
  int int_tmp;
  // declare the return object
  VISMatrix matrix;
  // get the number rows and cols
  good=true;
  if (VPF::set(r, pf[name][0]) != VPF::VALID)
    good=false;
  if (VPF::set(c, pf[name][1]) != VPF::VALID)
    good=false;
  if (!good)
    {
      cout << "ImplicitModel::bad matrix " << name << endl;
      return VISMatrix();
    }
  // create the return object
  matrix = VISMatrix(r,c);
  // read the data
  for (j = 0; j < r*c; j++)
    if (VPF::set(float_tmp, pf[name][j+2]) == VPF::VALID)
      matrix.poke(j/c, j%c) = float_tmp;
    else 
      if (VPF::set(int_tmp, pf[name][j+2]) == VPF::VALID)
	matrix.poke(j/c, j%c) = (float)int_tmp;
      else
	good=false;

  if (!good)
    {
      cout << "ImplicitModel::bad matrix " << name << endl;
      return VISMatrix();
    }
  return(matrix);
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
    if (VPF::set(float_tmp, pf[name][j+1]) == VPF::VALID)
      vector.poke(j) = float_tmp;
    else 
      if (VPF::set(int_tmp, pf[name][j+1]) == VPF::VALID)
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


main (int argc, char **argv)
{
  VISVolume<float> test, ref;
  VISVolumeFile vol_file;
  VISIm im;
  
  //declare VPF
  VPF::ParameterFile pf(argv[1]);

#ifdef DEBUG
  cout << "Reading in the parameter file" << endl << endl;
#endif
  
  /* set the registration method */
  

  // set the metric
  MutualInformation<AffineTransform> *metric = \
    new MutualInformation<AffineTransform> ();

  // set the Levenberg-Marquardt registrator
  Registrator< MutualInformation<AffineTransform>, \
               AffineTransform, Marquardt<MutualInformation<AffineTransform> > >\
    *regLM  = new Registrator< MutualInformation<AffineTransform>, \
    AffineTransform, Marquardt<MutualInformation<AffineTransform> > > ();

  //set the Levenberg-Marquardt optimizer
  Marquardt<MutualInformation<AffineTransform> > *optimizerLM \
    = new Marquardt<MutualInformation<AffineTransform> >();
  regLM->setOptimizer(optimizerLM);


  //set the Gradient Descent registrator
  Registrator< MutualInformation<AffineTransform>, \
               AffineTransform, GradientDescent<MutualInformation<AffineTransform> > >\
	       *regGD  = new Registrator< MutualInformation<AffineTransform>, \
    AffineTransform, GradientDescent<MutualInformation<AffineTransform> > > ();
  
  //set the Gradient Descent optimizer
  GradientDescent<MutualInformation<AffineTransform> > *optimizerGD \
    = new GradientDescent<MutualInformation<AffineTransform> >();
  regGD->setOptimizer(optimizerGD);  


  //set the transform
  AffineTransform *transform = new AffineTransform();

  int nA = 50;
  //set the number of samples for A and B
  if (VPF::set(nA, pf["SAMPLES"][0]) == VPF::VALID)
  {
    metric->setNA(nA);
    metric->setNB(nA);
  }
  else
    cout << "Warning! invalid number of samples" << endl;
  //set the metric
  regLM->setMetric(metric);
  regGD->setMetric(metric);

  //read in the reference volume
  string refName;
  if (VPF::set(refName, pf["REF_FILE"][0]) == VPF::VALID)
  {
    cout << "reading in the reference volume:" << endl << refName.c_str() << endl;
    ref = vol_file.read_float(refName.c_str());
    cout << "read in reference volume." << endl;
    regLM->setReferenceVolume(ref);
    regGD->setReferenceVolume(ref);
  }
  else
    cout << "Warning! invalid reference image filename specified" << endl;

  //set the number of images to be registered
  int numRegs;
  if (!(VPF::set(numRegs, pf["NUM_REGS"][0]) == VPF::VALID))
    cout << "Warning! the number of files to be registered was not specified" << endl;

  //read in the output directory
  string saveDir;
  if (!(VPF::set(saveDir, pf["OUTPUT_DIR"][0]) == VPF::VALID))
    cout << "Warning! improperly specified output directory" << endl;

  //set the number of subsampling levels to go through
  int numSubSamples=1;
  if (!(VPF::set(numSubSamples, pf["SUBSAMPLES"][0]) == VPF::VALID))
    cout << "Warning! did not specify the number of subsampled images" //
      " defaulting to one" << endl; 

  //read in the lambda value for levenberg-marqardt
  float lambdaLM = 0.0001;
  if (VPF::set(lambdaLM, pf["LM_LAMBDA"][0]) == VPF::VALID)
    regLM->setStepSize(lambdaLM);
  else
    cout << "Warning! invalid learning rate specified for Levenberg-Marquardt" << endl;
  
  int runs = 100;
  if (!(VPF::set(runs, pf["LM_RUNS"][0]) == VPF::VALID))
    cout << "Warning! invalid number of runs specified for Levenberg-Marquardt." <<
      endl << "Defaulting to one run" << endl;

  //read in the lambda value for gradient descent
  float lambdaGD = 0.01;
  if (VPF::set(lambdaGD, pf["GD_LAMBDA"][0]) == VPF::VALID)
    regGD->setStepSize(lambdaGD);
  else
    cout << "Warning! invalid learning rate specified for Gradient Descent" << endl;

  //by default set LM to large iters, it will stop on its own
  int iters=5000;
  //read in the number of iterations for gradient descent
  if (VPF::set(iters, pf["LM_ITERS"][0]) == VPF::VALID)
    regLM->setIters(iters);
  else
    cout << "Warning! invalid number of iterations specified" << endl;
  regLM->setIters(iters);
  //read in the number of iterations for gradient descent
  if (VPF::set(iters, pf["GD_ITERS"][0]) == VPF::VALID)
    regGD->setIters(iters);
  else
    cout << "Warning! invalid number of iterations specified" << endl;

  float Phi_v = 0.4;
  //read in the 1D covariance
  if (VPF::set(Phi_v, pf["PHI_V"][0]) == VPF::VALID)
    metric->setPhi_v(Phi_v);
  else
    cout << "Warning! invalid 1D covariance specified" << endl;

  boolean good = true, alsoGood = true;
  
  VISMatrix Phi_uv;
  //read in the 2D covariance 
  Phi_uv = readMatrix(pf, "PHI_UV", good); 
  if (good)
    metric->setPhi_uv(Phi_uv);
  else
    cout << "Warning! invalid covariance for u,v specified" << endl;

  //read in the initial transform
  VISVector params = readVector(pf, "INIT_TRANS", good);
  if (good)
    transform->setParams(params);
  else
    cout << "Warning! invalid transform specified" << endl;
  regLM->setInitialTransform(transform);
  
  // read in the weights that determine the amount of movement for the 
  // optimizer
  VISVector paramWeights = readVector(pf, "PARAM_WEIGHTS", good);
  if (good)
    metric->setParamWeights(paramWeights);
  else
    cout << "Warning! invalid parameter weights specified" << endl;
  
  //Read in the lower point for the region of interest
  VISVector lower = readVector(pf, "LOWER_ROI", good);
  //Read in the upper point for the region of interest
  VISVector upper = readVector(pf, "UPPER_ROI", alsoGood);
  if (good && alsoGood)
  {
    regLM->setROI(lower, upper);
    regGD->setROI(lower, upper);
  }
  else 
  {
    cout << "Warning! ROI incorrectly specified, defaulting"
	 << " to the entire image" << endl;
    lower = VISVector(2); lower = 0.0f;
    upper = VISVector(2); 
    upper.poke(0) = ref.width(); upper.poke(1) = ref.height();
  }
    
#ifdef DEBUG
  cout << "Finished reading in the parameter file, beginning registration" << endl << endl;
#endif  
 
  //begin the registration process


  //register for each test file
  char whichTest[350];
  string testName;
  char output[800];
  string viewerFile;
  for (int slice=0; slice<numRegs; slice++)
  {

    //read in the current test filename
    sprintf(whichTest, "TEST_FILE_%d", slice);
    if (VPF::set(testName, pf[whichTest][0]) == VPF::VALID)
    {
      cout << "reading in the test volume, " << testName << "." << endl;
      test = vol_file.read_float(testName.c_str());
      cout << "read in test volume." << endl << endl;
      regLM->setTestVolume(test); 
      regGD->setTestVolume(test); 
    }
    else
      cout << "Warning! invalid test image filename specified" << endl;

    //reset some of the registration parameters
    ///////////////////////////////////////////

    //reset the initial transform
    VISVector params = readVector(pf, "INIT_TRANS", good);
    if (good)
      transform->setParams(params);
    else
      cout << "Warning! invalid transform specified" << endl;
    regLM->setInitialTransform(transform);


    //reset the lambda value of the Levenberg-Marquardt Optimizer
    regLM->setStepSize(lambdaLM);
    optimizerLM->restart();
    
    //reset the Gradient Descent Optimizer
    optimizerGD->restart();
        
#ifdef DEBUG
    cout << "begining the Levenberg-Marquardt registration process" << endl;
#endif

    //begin the registration process
    //first register with Levenberg-Marquardt Optimization
    for (int k=0; k<runs; k++)
    {
      transform = regLM->registrate();
      cout << "transform: " << endl;
      transform->print();
      transform = regLM->getBestTransform();
      cout << "best transform: " << endl;
      transform->print();
    }

#ifdef DEBUG 
    cout << "begining the Gradient Descent optimization" << endl;
#endif
  
    //then register with Gradient Descent Optimization
    regGD->setInitialTransform(transform);
    transform = regGD->registrate();
    cout << "transform" << endl;
    transform->print();
    transform = regGD->getBestTransform();
    cout << "best transform: " << endl;
    transform->print();
 
    //make parameter file to send to the viewer
    cout << "writing viewing file for registration number " << slice << endl;
    ofstream viewer;
    //create viewerFilename
    sprintf(output, "%s/viewerParams%d.txt", saveDir.c_str(), slice);
    viewer.open(output);
    cout << "opened viewer: " << endl << output << endl << endl;
    viewer << "(REF_FILE \"" << refName << "\")" << endl << endl;
    viewer << "(TEST_FILE \"" << testName << "\")" << endl << endl;
    viewer << "(TRANS " << (regGD->getBestTransform())->getParams() << endl << ")" << endl << endl;
    viewer << "(LOWER_ROI " << lower << ")" << endl << endl;
    viewer << "(UPPER_ROI " << upper << ")" << endl << endl;
    viewer << "(CHECKS " << 17 << ")" << endl << endl;
    viewer << "(DOTS " << 7 << ")" << endl << endl;
    viewer << "(OUTPUT_DIR \"" << saveDir.c_str() << "\")" << endl << endl;
    viewer << "// metric = " << regGD->getBestMetric() << endl;
    viewer.close();
    
  }
  delete optimizerLM;
  delete optimizerGD;
  delete metric;
  delete transform;
  delete regLM;
  delete regGD;
  
  cout << "finished registering the images via MIregister" << endl;
  
}
 
   


