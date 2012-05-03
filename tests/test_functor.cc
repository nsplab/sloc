#include <iostream>

using namespace std;

// ----------------------------------------------------------------------------

class AdditiveFunctor
{
public:
    AdditiveFunctor(int x) : x(x) {}
    int operator()(int y) { return x + y; }
private:
    int x;
};

// ----------------------------------------------------------------------------

int main()
{
    AdditiveFunctor increment(1);
    cout << "increment(10) = " << increment(10) << endl;
    return 0;
}

