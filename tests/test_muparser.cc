// http://muparser.sourceforge.net/mup_version.html#idExample
#include <iostream>
#include <muParser.h>

using namespace std;

// function callback
double sqr(double a)
{
    return a * a;
}

int main(int argc, const char *argv[])
{
    try 
    {
        double a_val = 1;
        mu::Parser p;
        p.DefineVar("a", &a_val);
        p.DefineFun("sqr", sqr);
        p.SetExpr("sqr(a) * _pi + min(10,a)");

        for (std::size_t a = 0;  a < 10; ++a)
        {
            a_val = a;
            cout << p.Eval() << endl;
        }
    }
    catch (mu::Parser::exception_type& e)
    {
        cerr << e.GetMsg() << endl;
    }

    return 0;
}
