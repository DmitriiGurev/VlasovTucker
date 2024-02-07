#include "full.h"

#include <iostream>

using namespace VlasovTucker;
using namespace std;

int main()
{
    Tensor3d tensor(3, 3, 3);
    tensor.setConstant(1);

    Full full(tensor);

    cout << full << "\n";
    cout << full.Sum() << "\n";

    cout << (full + full) * (full + full) << "\n";
}