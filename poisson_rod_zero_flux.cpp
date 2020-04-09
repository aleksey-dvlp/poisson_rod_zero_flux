/*
u'' + u' = x*exp(x)
u(0) = 5
du(1)/dx = 0

Analitical solution: u(x) = 0.25*(exp(x)*(2*x-3) + exp(2-x)-exp(2)+23)

n = 1:
c[2] -2*c[1]+c[0]      c[2] - c[0]
-----------------   + ------------- = h*(exp(h)),  c[0] = 5
      h*h                 2*h

n = 2, N-1

c[n+1] -2*c[n]+c[n-1]      c[n+1] - c[n-1]
---------------------  + ------------------ = n*h*(exp(n*h))
      h*h                         2*h

n = N
 2*(c[N-1]-c[N])    
----------------  = N*h*(exp(N*h))
      h*h        
means that the second boundary condition (du(1)/dx=0) approximized by
adding a fictive node at the coordinate Nh+h and putting there a value
c[N+1] = c[N-1]
 */

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
using namespace std;


class ClassPoissRod {
  static constexpr int N = 70;
  const double L = 1.;
  const double h = L/N;
  double u[N+1] = {0.};
  double uan[N+1] = {0.}; //analitical solution
  ofstream file;

public:
  ClassPoissRod();
  ~ClassPoissRod();
  void norm();
  void solve();
  void analitsol();
  void writefile();
};

ClassPoissRod::ClassPoissRod() {
  auto data = "data";
  auto nstring = to_string(N);
  auto outstring = ".out";
  auto filename = data+nstring+outstring;
  file.open(filename);
}

ClassPoissRod::~ClassPoissRod() {
  file.close();
}

//analitical solution
void ClassPoissRod::analitsol() {
  for (int i = 0; i<=N; i++) {
    const double x = i*h;
    uan[i] = 0.25*(exp(x)*(2.*x-3.)+exp(2.-x)-exp(2.)+23.);
  }
}

void ClassPoissRod::writefile() {
  for (int i=0; i<=N; i++) {
    file << setprecision(10)<<i*h<<' '<<u[i]<<' '<<uan[i]<<
      ' '<<endl;
  }
}

//infinity norm
void ClassPoissRod::norm() {
  double s = 0.0;
  for (int n = 0; n<=N; n++) s = max(s, abs(u[n] - uan[n]));
  cout << "the infinity norm is "<<scientific<<setprecision(10)<<s<<endl;

}

  
  
void ClassPoissRod::solve() {
  double a[N+1], b[N+1], c[N+1], d[N+1], alpha[N+1], beta[N+1];
  const double invh2 = 1./h/h;
  const double invh  = 1./h;
  u[0] = 5.; //left Dirichlet boundaty condition
  a[1] = 0.0;
  b[1] = -2.*invh2;
  c[1] = 1.*invh2 + 0.5*invh;
  d[1] = h*exp(h)-u[0]*invh2 + 0.5*u[0]*invh;
  for (int n=2; n<N; n++) {
    a[n] = 1.*invh2 - 0.5*invh;
    b[n] = -2.*invh2;
    c[n] = 1.*invh2 + 0.5*invh;
    d[n] = n*h*exp(n*h);
  }
  a[N] = 2.*invh2;
  b[N] = -2.*invh2;
  c[N] = 0.;
  d[N] = N*h*exp(N*h);

  for (int n =1; n<=N; n++) {
    alpha[n] = c[n]/(b[n] - a[n]*alpha[n-1]);
    beta[n] = d[n] - a[n]*beta[n-1];
    beta[n] /= b[n] - a[n]*alpha[n-1];
  }
  for (int n=N; n>0; n--) u[n] = beta[n] - alpha[n]*u[n+1];
}

int main() {

  ClassPoissRod poissrod;
  poissrod.solve();
  poissrod.analitsol();
  poissrod.writefile();
  poissrod.norm();
  

}
  
