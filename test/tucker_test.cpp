#include "tucker.h"
#include <iostream>
#include <chrono>

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

using namespace VlasovTucker;
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
    // Matrix<double, 3, 9, RowMajor> A0;
    // A0 << 0.9073, 0.7158, -0.3698, 1.7842, 1.6970, 0.0151, 2.1236, -0.0740, 1.4429,
    //   0.8924, -0.4898, 2.4288, 1.7753, -1.5077, 4.0337, -0.6631, 1.9103, -1.7495,
    //   2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335, -0.2716;

    // Tensor<double, 3> tensor = Folding3D(A0, 0, 3, 3, 3);

    // Tucker tucker(tensor);

    // cout << tucker.Reconstructed() << "\n\n";

    // VectorXd u(3);
    // u << 2, 2, 2;

    // cout << u << "\n\n";

    // Tucker t2({u, u, u});

    // // cout << (tucker / t2).Reconstructed() << "\n\n";

    // cout << tucker.Reconstructed() << "\n\n";

    //cout << Reflection(tucker, 1).Reconstructed() << "\n\n";

    //cout << (tucker * tucker).Reconstructed() << "\n\n";

    //cout << tucker.Reconstructed() << "\n\n";

    //Tucker tucker1(tensor);
    //std::cout << tucker1 << "\n\n";
    //Tucker tucker2(tensor);
    //tucker2.Compress(0.1, 2);
    //std::cout << tucker2.Reconstructed() << "\n\n";

    //Tucker tucker(tensor);
    //std::cout << A0.norm() << " " << tucker1.norm();

    //cout << tensor(0,1,2) << " " << (tucker*2-tucker).At(0, 1, 2);

    // int n = 3;
    // Tensor<double, 3> tensorLarge(n, n, n);
    // tensorLarge.setRandom();
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n; j++) {
    //         for (int k = 0; k < n; k++) {
    //             double x = i / (double)n; 
    //             double y = j / (double)n; 
    //             double z = k / (double)n; 

    //             double r = sqrt(pow(0.4 - x, 2) + pow(0.3 - y, 2) + pow(0.5 - z, 2));
                
    //             // tensorLarge(i, j, k) = exp(-pow(r / 0.2, 2)); 
    //             tensorLarge(i, j, k) = abs(x - 0.5) + abs(2 * y - 0.4) + 1.5 * z - 0.6 + z * z; 
    //         }
    //     }
    // }

    // // cout << tensorLarge << "\n";

    // Tucker tuckerLarge(tensorLarge, 1.0e-5);
    // cout << tuckerLarge.Core() << "\n\n";
    // for (auto u : tuckerLarge.U())
    //     cout << u << "\n\n";
    // cout << tuckerLarge << "\n\n";

    // cout << (tuckerLarge.Reconstructed() - tensorLarge).square().sum().sqrt() << "\n";

    // Tensor<double, 3> tensor(10, 10, 10);
    // tensor.setRandom();

    // Tucker r1(tensor, 0, 1);
    // Tucker r(tensor, 0, 5);

    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // for (int i = 0; i < 100; i++)
    // {
    //     Tucker term1 = r + r;
    //     term1.Compress(0, 5);

    //     Tucker term2 = r * term1;
    //     term2.Compress(0, 5);

    //     Tucker term3 = term2 + term2;
    // }

    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1.0e6 << " [s]" << std::endl;

    // for (int i = 0; i < 10000; i++)
    // {
    //     Tucker tucker(tensor);
    //     tucker.Compress(1e-6, 2);
    //     // tucker = Tucker(tucker.Reconstructed(), 1e-6, 2);
    // }

    // for (int maxR : {4, 3, 2, 1})
    // {
    //     tucker.Compress(0, maxR);
    //     cout << tucker << "\n";
    //     cout << (tucker.Reconstructed() - tensor).square().sum().sqrt() << "\n\n";
    // }

    // for (double eps : {1e-10, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 0.8, 1.0, 1.5})
    // {
    //     tucker.Compress(eps);
    //     cout << tucker << "\n";
    //     cout << (tucker.Reconstructed() - tensor).square().sum().sqrt() << "\n\n";
    // }

    Tensor<double, 3> tensor(5, 3, 3);
    tensor.setRandom();

    Tucker tucker(tensor, 1e-6, 4);
    // for (int i = 0; i < 100; i++)
    // {
    //     cout << "i = " << i << "\n";

    //     Tucker t1 = tucker;
    //     tucker += t1;

    //     Tucker t2 = 0.5 * tucker;
    //     tucker -= t2;

    //     tucker.Compress(1e-6);

    //     cout << tucker << "\n";
    //     cout << tucker.Reconstructed() << "\n\n";
    // }

    for (int i = 0; i < 10; i++)
    {
        tucker += tucker;
        tucker -= 0.5 * tucker;
        tucker += tucker;
        tucker -= 0.5 * tucker;
        // cout << tucker << "\n";
        tucker.Compress(1e-6, 4);
        // cout << tucker << "\n";
        cout << tucker.Reconstructed() << "\n\n";
    }
}