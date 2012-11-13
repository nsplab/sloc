#include <sloc/hptimer.h>
#include <unistd.h>
#include <iostream>

using namespace std;

inline void ms_sleep(int ms)
{
    usleep(1000 * ms);
}

int main(void)
{
    unsigned long N = 1000;
    sloc::hptimer* T = sloc::hp::get_timer("One-thousand calls to ms_sleep(100)");
    
    cout << T->get_name() << endl;
    for (unsigned long i = 1; i <= N; i++)
    {
        T->start();
        cout << "  " << i << "\r" << std::flush;
        ms_sleep(100);
        T->stop();
    }
    cout << "    \r" << std::flush;

    sloc::hp::report();

    return 0;
}
