#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

#define Nx 50
#define Ny 50
#define repeats 30
#define steps 500000
#define Tmin 0.5
#define Tmax 5.0
#define Jval 1.0
typedef float real;

inline size_t bc(size_t i, size_t N) {
  // Apply the boundary condition:
  if (i == -1) {
    return N - 1;
  }
  else if (i == N) {
    return 0;
  }
  else {
    return i;
  }
}

inline real dE(size_t i, size_t j, real m[Nx][Ny], real J, real T, real H) {
  return - 2 * J * m[i][j] * (
              m[bc(i-1, Nx)][j] +
              m[bc(i+1, Nx)][j] +
              m[i][bc(j-1, Ny)] +
              m[i][bc(j+1, Ny)]
            ) -
            m[i][j] * H;
}

inline void step(real m[Nx][Ny], real T, real J, real H) {
  size_t sx = rand() % Nx;
  size_t sy = rand() % Ny;
  real deltaE = - dE(sx, sy, m, J, T, H);
  if (deltaE < 0) {
    m[sx][sy] *= -1;
  }
  else if (exp(-deltaE / T) > ((float) rand())/RAND_MAX) {
    m[sx][sy] *= -1;
  }
}

void write_field(real m[Nx][Ny], size_t run, double T, double J) {
  std::string fname = "run_" + std::to_string(run) + "_" + std::to_string(Nx) + "_" + std::to_string(Ny);
  fname += "_T_" + std::to_string(T) + "_J_" + std::to_string(J) + ".txt";
  std::ofstream fout(fname);

  for(size_t i = 0; i < Nx; i++) {
    for(size_t j = 0; j < Ny; j++) {
      fout << m[i][j] << ",";
    }
    fout << std::endl;
  }
}

int main() {
  srand(time(NULL));
  real m[Nx][Ny];
  real J = 1.0;
  real H = 0.01;
  // Initialise

  for(int r = 0; r < repeats; r++) {
    std::cout << "Run: " << r << std::endl;
    for(int i = 0; i < Nx; i++) {
      for(int j = 0; j < Ny; j++) {
        m[i][j] = 2*(rand() % 2) - 1;
      }
    }

    for (double T = 8; T > 0.01; T -= 0.05) {
      for(int s = 0; s < steps; s++) {
        step(m, T, J, H);
      }
      write_field(m, r, T, J);
    }
  }



  return 0;







}
