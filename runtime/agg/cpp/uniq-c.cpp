#include "main.h"
#include <map>

void aggregate() noexcept
{
    size_t padding = 0;
    if (input1())
    {
        size_t count;
        input1() >> count;
        padding = input1().tellg();
        input1().seekg(std::ios_base::beg);
    }

    std::string commonKeyCand{};
    size_t delta = 0;
    bool outputed = false;
    if (input2())
    {
        input2() >> delta;
        padding = input2().tellg();
        std::getline(input2(), commonKeyCand);
    }

    {
        std::string key;
        size_t count;
        while (true)
        {
            input1() >> count;
            if (!input1())
                break;
            std::getline(input1(), key);
            if (key == commonKeyCand)
            {
                count += delta;
                outputed = true;
            }
            
            output().width(padding);
            output() << count << key << '\n';
        }
    }
    if (!outputed)
    {
        output().width(padding);
        output() << delta << commonKeyCand << '\n';
    }

    output() << input2().rdbuf();
}
