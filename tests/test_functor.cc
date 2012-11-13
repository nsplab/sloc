#include <iostream>
#include <sloc/utils.h>

class AdditiveFunctor 
{
public:
    AdditiveFunctor(int x) : x(x) {}
    int operator()(int y) { return x + y; }
private:
    int x;
};

int main(void)
{
    AdditiveFunctor increment(1);
    _PRINT_VALUE(increment(10)); 
    return 0;
}
