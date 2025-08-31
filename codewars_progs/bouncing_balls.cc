#include "all_codewars.h"

using namespace std;
int codewars_mysol::bouncingBall(double h, double bounce, double window)
{
  if ((bounce > 0.0 && bounce < 1.0) && (window < h) && h > 0.0)
  {
    int count = 1;
    double h1 = h;
    while (h1 > window)
    {
       h1 *= bounce;
       count += 2;
    }
    count -= 2;
    return count;
  } else
  {
    return -1;
  }
}
