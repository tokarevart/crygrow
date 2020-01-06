#pragma once
#include <cstdlib>
#include <cmath>


namespace itr {

enum class direction {
    forward,
    reverse
};

// [begin; end)
// if forward : begin, begin+1, ... end-1
// else       : end-1, end-2, ... begin
template <typename Signed>
class iteration {
public:
    Signed begin() const {
        return m_start;
    }
    Signed end() const {
        return m_endm1 + 1;
    }
    Signed next() {
        return m_begin + std::abs(m_start - ++m_current);
    }
    bool has_ended() const {
        return m_current > m_endm1;
    }

    void init(Signed begin, Signed end, direction dir) {
        m_start = dir == direction::forward ? begin : end - 1;
        m_begin = begin;
        m_current = begin;
        m_endm1 = end - 1;
    }

    iteration(Signed begin, Signed end, direction dir = direction::forward) {
        init(begin, end, dir);
    }

private:
    Signed m_begin;
    Signed m_current;
    Signed m_endm1;
    Signed m_start;
};

} // namespace itr
