#include <iostream>
#include <memory>

bool stream_copy_n(std::istream& in, std::ostream& out, std::size_t count) noexcept
{
    const std::size_t buffer_size = 256 * 1024;
    std::unique_ptr<char[]> buffer = std::make_unique<char[]>(buffer_size);
    while(count > buffer_size)
    {
        in.read(buffer.get(), buffer_size);
        out.write(buffer.get(), buffer_size);
        count -= buffer_size;
    }

    in.read(buffer.get(), count);
    out.write(buffer.get(), count);

    return in.good() && out.good();
}