#include <sloc/progress_timer.h>
#include <unistd.h>

int main(int argc, const char *argv[])
{
    int N = 10000;
    sloc::ProgressTimer ptimer;

    std::cerr << ptimer.header("cycles");
    ptimer.start(N);
    for (int i = 1; i <= N; i++)
    {
        usleep(500);
        std::cerr << ptimer.update(i);
    }
    ptimer.update(N);
    std::cerr << std::endl;

    return 0;
}
