#include "base.h"
#include "util.h"
#include "vec.h"
#include "mat.h"
#include "optim.h"
#include "redund.h"

using namespace std;
int main()
{
	mat<double> A(3,3);
	mat<double> B(3,3);

	for (std::size_t Y =0; Y<A.size_y;Y++){
		for (std::size_t X =0; X<A.size_x;X++){
			A(Y,X) = 1 + Y + X*A.size_y;
		}
	}

	for (std::size_t Y =0; Y<B.size_y;Y++){
		for (std::size_t X =0; X<B.size_x;X++){
			B(Y,X) = 7 + Y + X*B.size_y;
		}
	}

	mat<double> C(2,2);
	C.set_all(10);
	C.print();
	mat<double> D;
	D=C;
	D.print();
	mat<bool>E = D;
	E.print();
	mat<int> G(E);
	G.print();
	mat<int> H = G;
	H.print();
	cout << "OPERATORS" << endl;
	cout << "+ && -" << endl;
	(-C).print();
	(+C).print();
	cout << "+ && -T" << endl;
	(A+10.0).print();
	(B-123.0).print();
	cout << "+ && -mat<T>" << endl;
	(A+B).print();
	(A-B).print();
	cout << "%" << endl;
	(A%B).print();
	cout << "+= && -= T" << endl;
	C.set_all(12);
	C+=3;
	C.print();
	C-=10;
	C.print();
	cout << "+= && -= mat<T>" << endl;
	C+=D;
	C.print();
	C-=D;
	C.print();
	cout << "*= && /= T" << endl;
	mat<double> v15(3,2);
	C*=3;
	C.print();
	C/=10;
	C.print();
	cout << "*= && /= mat<T>" << endl;
	C*=C;
	C.print();
	C/=C;
	C.print();
	cout << "det && adj && inv" << endl;
	mat<double> M(2,2);
	M(0,0,4);	M(0,1,7);
	M(1,0,2);	M(1,1,6);
	cout << "det" << endl;
	cout << M.det() << endl;
	cout << "adj" << endl;
	(M.adjoint()).print();
	cout << "inv" << endl;
	(M.inverse()).print();
	cout << "concat" << endl;
	(M.concat(M)).print();
	cout << "row" << endl;
	M.row(0).print();
	cout << "col" << endl;
	M.col(0).print();
	cout << "trunc" << endl;
	M.print();
	M.trunc(0,0,2,1).print();
	cout<< "load mat from file" << endl;
	mat<double> vf;
	mat<double> N(5,3);
	N(0,0, 0);N(0,1, 1);N(0,2, 0);
	N(1,0, 0);N(1,1, 0);N(1,2, 1);
	N(2,0, 1);N(2,1,-1);N(2,2, 0);
	N(3,0, 1);N(3,1, 0);N(3,2,-1);
	N(4,0, 1);N(4,1, 0);N(4,2,-1);
	mat<double> L =N.redund();
	cout << "redunding:" << endl;
	L.print();
	vec<double> obj(L.size_x);
	obj.set_all(1); obj(0) = 0;
	vec<double> sol = L.lpsolve(obj);
	cout << "optim:" << endl;
	sol.print();
	mat<double> minimal("temp2.ine");
	cout << minimal.size_y << "," << minimal.size_x << endl;
	minimal.print();
	mat<double> Q = A % B;
	cout << "A" << endl;
	A.print();
	cout << "B" << endl;
	B.print();
	cout << "C = A % B" << endl;
	Q.print();
	mat<double> ve("eig.test");
	mat<double> eig = ve.jacobi_eigen();
	eig.print();

	mat<double> load_from_file = from_ine<double>("000.ine");
	load_from_file.print();
}
