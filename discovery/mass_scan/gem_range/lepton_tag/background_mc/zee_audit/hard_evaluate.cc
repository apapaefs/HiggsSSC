#include "CPPProcess.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    CPPProcess process;
    process.initProc(argv[1]);
    double p[5][4];
    std::vector<double*> momenta{p[0],p[1],p[2],p[3],p[4]};
    while (std::cin >> p[0][0]) {
        for (int particle = 0; particle < 5; ++particle)
            for (int k = (particle == 0); k < 4; ++k)
                if (!(std::cin >> p[particle][k])) return 2;
        process.setMomenta(momenta);
        process.sigmaKin();
        std::cout << "M2 " << std::setprecision(17) << process.getMatrixElements()[0] << '\n';
    }
    return 0;
}
