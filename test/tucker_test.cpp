#include "tucker.h"
#include <iostream>

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

using namespace Eigen;
using std::cout;

#include <iostream>

Tensor<double, 3> Folding3D(MatrixXd unfolding, int index, int I0, int I1, int I2)
{
    Tensor<double, 3> folding(I0, I1, I2);

    switch (index)
    {
    case 0:
        for (int i0 = 0; i0 < I0; i0++)
        {
            for (int i1 = 0; i1 < I1; i1++)
            {
                for (int i2 = 0; i2 < I2; i2++)
                {
                    folding(i0, i1, i2) = unfolding(i0, i2 + i1 * I2);
                }
            }
        }
        break;
    case 1:
        for (int i1 = 0; i1 < I1; i1++)
        {
            for (int i2 = 0; i2 < I2; i2++)
            {
                for (int i0 = 0; i0 < I0; i0++)
                {
                    folding(i0, i1, i2) = unfolding(i1, i0 + i2 * I0);
                }
            }
        }
        break;
    case 2:
        for (int i2 = 0; i2 < I2; i2++)
        {
            for (int i0 = 0; i0 < I0; i0++)
            {
                for (int i1 = 0; i1 < I1; i1++)
                {
                    folding(i0, i1, i2) = unfolding(i2, i1 + i0 * I1);
                }
            }
        }
        break;
    }
    return folding;
}

int main()
{
    Matrix<double, 3, 9, RowMajor> A0;
    A0 << 0.9073, 0.7158, -0.3698, 1.7842, 1.6970, 0.0151, 2.1236, -0.0740, 1.4429,
      0.8924, -0.4898, 2.4288, 1.7753, -1.5077, 4.0337, -0.6631, 1.9103, -1.7495,
      2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335, -0.2716;

    Tensor<double, 3> tensor = Folding3D(A0, 0, 3, 3, 3);

    Tucker tucker(tensor);

    cout << tucker.Reconstructed() << "\n\n";

    VectorXd u(3);
    u << 2, 2, 2;

    cout << u << "\n\n";

    Tucker t2({u, u, u});

    cout << (tucker / t2).Reconstructed() << "\n\n";

    //cout << tucker.Reconstructed() << "\n\n";

    //cout << Reflection(tucker, 1).Reconstructed() << "\n\n";

    //cout << (tucker * tucker).Reconstructed() << "\n\n";

    //cout << tucker.Reconstructed() << "\n\n";

    //Tucker tucker1(tensor);
    //std::cout << tucker1 << "\n\n";
    //Tucker tucker2(tensor);
    //tucker2.Recompress(0.1, 2);
    //std::cout << tucker2.Reconstructed() << "\n\n";

    //Tucker tucker(tensor);
    //std::cout << A0.norm() << " " << tucker1.norm();

    //cout << tensor(0,1,2) << " " << (tucker*2-tucker).At(0, 1, 2);
}