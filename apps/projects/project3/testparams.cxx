#include <array.h>
#include <param.h>
#include <imagefile.h>
#include <matrix.h>
#include <image.h>



//
// A helper function which reads, from a ParameterFile object, a
// matrix of the form (name rows cols ....data in row major order...)
// Returns an invalid matrix if it fails (prints something)
//
VISMatrix readMatrix(VPF::ParameterFile pf, const char* name)
{
  int r, c, j;
  float float_tmp;
  int int_tmp;
  // declare the return object
  VISMatrix matrix;
  // get the number rows and cols
  boolean good=true;
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


//
// Program reads in a parameter file to specify a mosaicing operation
// in accordance with Project 2.  The format of the file is ...
//
//CORRESPONDENCE_MATRIX <int> <int> <floats in row major order> Matrix
//(square) --- rows and cols are image numbers, and upper half of this
//matrix gives correspondences available //IMAGE_FILE_<N> <string>
//Gives image file name for the Nth image
//
// CORRESPONDENCES_<M>_<N> <int> 2 <floats in row major order
// Specifies a Px2 matrix (p is # of correspondences) which give the
// control points in the Mth image that will match up with those in the 
// Nth image.  Note there must be a CORRESPONDENCES_<N>_<M> and it 
// must be the same size
// 
// TARGET_IMAGE <int> specify which target image coordinate system to
// produce the final result
main(int argc, char** argv)
{
  VISImageFile im_file;
  if (argc < 2)
    {
      cout << "usage: readparams parameterfile" << endl;
      exit(-1);
    }
  // create parameter file from the given argument
  VPF::ParameterFile pf(argv[1]);

  // read in the matrix of correspondences for which we will look  
  VISMatrix correspondence_matrix = readMatrix(pf, "CORRESPONDENCE_MATRIX");

  // print it out just for debugging
  cout << "correspondence matrix " << endl << correspondence_matrix << endl;

  // make sure it's valid and square
  if ((!correspondence_matrix.isValid())
      ||(!(correspondence_matrix.r() == correspondence_matrix.c())))
    {
      cout << "bad correspondence matrix in param file: " << argv[1] << endl;
      exit(-1);
    }

  // the number of images is the size of this matrix  (one row/col per
  // image)
  int num_images = correspondence_matrix.r();

  // create an 2D array that contains all of the control points
  VISArray< VISArray< VISMatrix > >  correspondences;

  int i, j, k;
  char keyword[80];
  VISMatrix pointsA, pointsB;


  // put empty matrices along the diagonal
  for (i = 0; i < num_images; i++)
    (correspondences.poke(i)).poke(i) = VISMatrix();

  // for all entries in the upper half
  for (i = 0; i < num_images; i++)
    for (j = i+1; j < num_images; j++)
      if (correspondence_matrix.peek(i, j))
      {
	// look for the two matching correpondence matrices
	sprintf(keyword, "%s_%d_%d", "CORRESPONDENCES", i, j);
	pointsA = readMatrix(pf, keyword);
	sprintf(keyword, "%s_%d_%d", "CORRESPONDENCES", j, i);
	pointsB = readMatrix(pf, keyword);
	// make sure they are the same size
	if (!pointsB.compareSize(pointsA))
	  {
	    cout << "bad control matrices " << i << " " << j << " in param file: " << argv[1] << endl;
	    exit(-1);
	  }
	// transpose so that individual points are columns
	// and put them in the 2D array of matrices
	(correspondences.poke(i)).poke(j) = pointsA.t();
	(correspondences.poke(j)).poke(i) = pointsB.t();
      }
      else
	{
	  // if those correspondences are not expectd, put
	  // in invalid matrices
	  (correspondences.poke(i)).poke(j) = VISMatrix();
	  (correspondences.poke(j)).poke(i) = VISMatrix();
	}

  int target_image;

  // get the target image number  
  if ((VPF::set(target_image, pf["TARGET_IMAGE"][0]) != VPF::VALID)
      ||(target_image >= num_images))
    {
      cout << "bad target image in param file: " << argv[1] << endl;
      exit(-1);
    }

  // print out the correspondeces --- for debugging
  for (i = 0; i < num_images; i++)
    for (j = i+1; j < num_images; j++)
      if (correspondence_matrix.peek(i, j))
      {
	cout << "The correspondences between image " << i 
	     << " and image " << j << " are: "  << endl;
	pointsA = correspondences[i][j];
	pointsB = correspondences[j][i];
	for (k = 0; k < pointsA.c(); k++)
	  {
	    cout << "c" << k << "=(" << pointsA.peek(0, k) << "," << pointsA.peek(1, k) << ")  : c'" << 
	      k << "=(" << pointsB.peek(0, k) << "," << pointsB.peek(1, k) << ") " << endl;
	  }
	cout << endl;
      }


  // read in the image and put them in an array
  VISArray< VISImage< float > > images;
  char filename[80];
    for (i = 0; i < num_images; i++)
      {
        sprintf(keyword, "%s_%d", "IMAGE_FILE", i);    
        if (VPF::set(filename, pf[keyword][0]) != VPF::VALID)
  	{
  	  cout << "no image file " << keyword << " in param file: " << argv[1] << endl;
  	  exit(-1);
  	}
        images.poke(i) = VISImage<float>(im_file.read(filename));
        if (!(images.peek(i)).isValid())
  	{
  	  cout << "bad image file " << filename << " in param file: " << argv[1] << endl;
  	  exit(-1);
  	}
      }

}
