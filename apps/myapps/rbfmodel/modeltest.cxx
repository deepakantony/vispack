#include "dipole.h"

main(int argc, char** argv)
{

  VISImageFile im_file;

  if (argc < 2)
    {
      cout << "usage: modeltest parameterfile" << endl;
      exit(-1);
    }

  VPF::ParameterFile pf(argv[1]);
  
  ImplicitModel model(pf);

  #define IMAGE_SIZE 100
  //  im_file.write(model.getImage(200,200), "modeltest_out.fit");
  im_file.write_jpeg(model.getImageModel(IMAGE_SIZE,IMAGE_SIZE), "outfile_0.jpg", 100);
  //  im_file.write(model.getImage(100,100).zeroCrossings(), 
  //  		"modeltest_zeros.fit");
  model.printCurvatures();
  //  im_file.write(model.getImageCurve(100,100), 
  //		"modeltest_curve.fit");
  int i;
  char filename[80];
  for (i = 1; i <= 40; i++)
    {
      model.deform();
      model.printCurvatures();
      if (i%1 == 0)
	{
	  sprintf(filename, "outfile_%d.jpg", i);
	  im_file.write_jpeg(model.getImageModel(IMAGE_SIZE, IMAGE_SIZE), 
			     filename, 100);
	  cout << "saved output " << filename << endl;
	}
    }
  //  im_file.write(model.getImage(100,100), "modeltest_out2.fit");

  //  im_file.write(VISImage<float>(model.getImageModel(100,100)), "model_test_color.fit");
  //  VISImageRGBA im_color = model.getImageModel(100,100);
  //VISImage<float> im_color = VISImage<float>(100,100);
  //  im_color.print();
  //  im_file.write_jpeg(im_color, "model_final.jpg", 100);

  //  model.printCurvatures();

  cout << "done" << endl;

}




