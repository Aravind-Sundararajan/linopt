#include "vec.h"
using namespace std;

double f(double n,double d){return pow(n,d);}
double g(double n){return f(n,.5);}
double h(double n){return f(n,2.0);}

int main()
{
	cout<< "load vec from file" << endl;
	vec<int> v1(10);
	v1.set_all(10);
	v1.print();
	vec<double> v2;
	v2=v1;
	v2.print();
	vec<double>v3 = v1;
	v3.print();
	vec<bool>v4;
	v4 = v1;
	v4.print();
	vec<int>v5(v2);
	v5.print();
	vec<int>v6 = v2;
	v5.print();
	cout << "OPERATORS" << endl;
	cout << "+ && -" << endl;
	vec<double> v7 = +v2;
	vec<double> v8 = -v2;
	v2.print();
	v7.print();
	v8.print();
	cout << "+ && -T" << endl;
	v2+10.0;
	vec<double> v9(v2+10.0);
	vec<double> v10(v2-123.0);
	v2.print();
	v9.print();
	v10.print();
	cout << "+ && -vec<T>" << endl;
	vec<double> v11(v10+v10);
	vec<double> v12(v11-v11);
	v11.print();
	v12.print();
	cout << "%" << endl;
	double sol = v11%v11;
	cout << sol << endl;
	cout << "+= && -= T" << endl;
	vec<double> v13(10);
	v13.set_all(12);
	v13+=3;
	v13.print();
	v13-=10;
	v13.print();
	cout << "+= && -= vec<T>" << endl;
	vec<double> v14(3);
	v14(0,1);
	v14(1,2);
	v14(2,3);
	v14+=v13;
	v14.print();
	v14-=v13;
	v14.print();
	cout << "*= && /= T" << endl;
	vec<double> v15(3);
	v15(0,1);
	v15(1,2);
	v15(2,3);
	v15*=3;
	v15.print();
	v15/=10;
	v15.print();
	cout << "*= && /= vec<T>" << endl;
	vec<double> v16(3);
	v16.set(0,1);
	v16.set(1,2);
	v16.set(2,3);
	v16*=v13;
	v16.print();
	v16/=v14;
	v16.print();
	cout << "test vecfun:" << endl;
	vec<double> v21(10);
	for (std::size_t X=0; X<v21.size;X++){
		v21(X,double(X));
	}
	v21.print();
	v21 = v21.fun(h);
	v21.print();
	v21 = v21.fun(g).fun(g).fun(g);
	v21.print();
	return 0;
};
