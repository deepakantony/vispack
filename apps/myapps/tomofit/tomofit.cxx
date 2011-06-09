
#include <matrix.h>
#include <projection.h>

#define size (40)
#define theta (M_PI/4.0f)
//#define theta (0.0f)
#define x0 (0.0f)
#define y0 (2.0f*size)
//#define scale (size/5.0f)
#define scale (30.0f)

main(int argc, char** argv)
{
  Projection proj(size);
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

  params.poke(0) = theta + 0.02;

}
