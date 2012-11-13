#ifndef HPTIMER_H
#define HPTIMER_H

#include <sys/time.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <sstream>
#include <map>

// ----------------------------------------------------------------------------

#define MAX_TIMEVAL (0x7fffffffffffffff)

typedef __int64_t hptime_t;

inline void tic(hptime_t& t)
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    t = tv.tv_sec * (hptime_t)1000000 + tv.tv_usec;
}

inline std::string format_hptime(double t)
{
    // XXX: equivalent of "%g" for output streams?
    std::ostringstream ss;
    if (t > 1e9)
        ss << (t * 1e9) << " ns";
    else if (t > 1e6)
        ss << (t * 1e6) << " us";
    else if (t > 1e3)
        ss << (t * 1e3) << " ms";
    else
        ss << t << " s";
    return ss.str();
}



namespace sloc {

// ----------------------------------------------------------------------------

class hptimer
{
public:
    hptimer(std::string name) : m_name(name)
    {
        reset();
    }

    ~hptimer()
    {
    }

    inline void start(void)
    {
        tic(m_start);
    }

    inline void stop(void)
    {
        tic(m_stop);
        m_cur = m_stop - m_start;
        if (m_cur < 0) return;
        m_total += m_cur;
        m_count++;
        if (m_min > m_cur) { m_min = m_cur; }
        if (m_max < m_cur) { m_max = m_cur; }
    }

    void reset(void)
    {
        m_start = m_stop = 0;
        m_cur = m_total = 0;
        m_min = MAX_TIMEVAL;
        m_max = -1;
        m_count = 0;
    }

    void report(void)
    {
        using namespace std;
        const double sec = 1e6;
        cout << "  Timer: " << m_name << " - called " << m_count << " times" << endl;
        if (m_count)
        {
            cout << "    Minimum: " << format_hptime(m_min / sec) << endl;
            cout << "    Maximum: " << format_hptime(m_max / sec) << endl;
            cout << "    Average: " << format_hptime(m_total / sec / m_count) << endl;
            cout << "    Minimum: " << format_hptime(m_min / sec) << endl;
        }
    }

    const std::string& get_name() { return m_name; }

private:
    hptime_t m_start, m_stop;
    hptime_t m_cur, m_min, m_max, m_total;
    int m_count;
    std::string m_name;
};

// ----------------------------------------------------------------------------

class hp
{
public:
    typedef std::map<std::string,hptimer*> hp_map_t;

    ~hp()
    {
        hp_map_t::iterator it;
        for (it = m_timers.begin(); it != m_timers.end(); ++it)
        {
            delete it->second;
            it->second = 0;
        }
        m_timers.clear();
    }

    static hptimer* get_timer(const char *name)
    {
        hp* instance = hp::get_instance();
        if (!instance->m_timers.count(name))
            instance->m_timers[name] = new hptimer(name);
        return instance->m_timers[name];
    }

    static hp* get_instance(void)
    {
        if (!m_instance)
            m_instance = new hp();
        return m_instance;
    }

    static void report(void)
    {
        hp* instance = hp::get_instance();
        hp_map_t::iterator it;
        for (it = instance->m_timers.begin(); it != instance->m_timers.end(); ++it)
        {
            hptimer* t = it->second;
            t->report();
        }
    }

private:
    // make singleton!
    hp() {}
    static hp* m_instance;

    // keep collection of timers
    hp_map_t m_timers;
};

// ----------------------------------------------------------------------------

} // namespace sloc

#endif
