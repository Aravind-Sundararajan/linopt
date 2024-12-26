#include "base.h"
#include "util.h"
#include "kvp.h"
#include "vec.h"
#include "mat.h"
#include "optim.h"
using namespace std;
int main()
{
	mat<double> H(2,3);
	H(0,0,1); H(0,1,-1); H(0,2,0);
	H(1,0,1); H(1,1,0); H(1,2,-1);
	H.print();
	vec<double> c(H.size_x);
	c(0,0); c(1,-1); c(2,-1);
	optim<double> o(H,c);
	vec<double>solution = (o.simplex()).trunc(0,H.size_x-1);
	solution.print();
	// c.set_all(0);
	// optim<double> o2(H,c);
	// vec<double>solution2 = o2.simplex();
	// solution2.print();
}
