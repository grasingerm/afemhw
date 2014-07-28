#include <armadillo>
#include <iostream>

static const int e_aux_data[12] =
{
    1, 0, 0,
    0, 0, 1,
    0, 0, 1,
    0, 1, 0
};

static const arma::Mat<int>::fixed<3,4> e(e_aux_data);

int main()
{
    std::cout << e << std::endl;
    return 0;
}
