#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
#include "kvp.h"
using namespace std;
int main()
{
	double f = 3;
	vec<double> v(3);
	v(0,1);v(1,10);v(2,20);
	kvp<double,vec<double>> p(f,v);
	v(0,11);v(1,101);v(2,201);
	cout << p.first << endl;
	p.second.print();
	kvp<vec<double>,double>  s =p.swap();
	cout << s.second << endl;
	s.first.print();
	v.print();
	kvp< kvp< vec<double> , double > , double > q(s,f);
	cout << q.second << endl;
	cout << q.first.second << endl;
	q.first.first.print();
}
