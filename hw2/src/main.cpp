#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
    double (*foo)(const double);

    foo = &yolo;

    cout << foo(3.4) << endl;

    return 0;
}

