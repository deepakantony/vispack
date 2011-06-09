
#include <matrix.h>
#include <projection.h>

#define size (40)
#define theta (M_PI/4.0f)
//#define theta (0.0f)
#define x0 (0.0f)
#define y0 (5.0f*size)
//#define scale (size/5.0f)
#define scale (10.0f)

main(int argc, char** argv)
{
  OrthoProjection proj(size);
  Square sq;
  VISVector params(4);
  params.poke(0) = theta;
  params.poke(1) = x0;
  params.poke(2) = y0;
  params.poke(3) = scale;

  sq.setParams(params);
  VISVector data = proj.project(sq);
  for (int i = 0; i < size; i++)
    {
      cout << data.peek(i) << " ";
    }
  cout << endl;

  cout << sq.getParams() << endl;
  params.poke(0) = theta - 0.05;
  params.poke(1) = x0 + 3.0f;
  params.poke(3) = 1.05*scale;
  //  params.poke(2) = y0 + 5.0f;
  sq.setParams(params);  
  SystemFitModel<OrthoProjection, Square> system(sq, proj, data);

  Solver solver;
  cout << system.getParams() << endl;
  solver.solve(system);
  cout << system.getParams() << endl;
  cout << "done" << endl;
}
