#pragma once
#include <cstdlib>
#include <cmath>


namespace itr {

enum class direction {
    forward,
    reverse
};

// begin < end:
// [begin; end)
//   if forward: begin, begin+1, ... end-1
//   else      : end-1, end-2,   ... begin
//
// begin > end:
// [end; begin)
//   if forward: begin, begin-1, ... end+1
//   else      : end+1, end+2,   ... begin
//
template <typename Signed>
class iteration {
public:
    Signed first() const {
        return m_begin;
    }
    Signed last() const {
        return m_end - m_step;
    }
    Signed begin() const {
        return m_begin;
    }
    Signed end() const {
        return m_end;
    }
    Signed next() {
        return m_current += m_step;
    }
    bool has_ended() const {
        return m_current == m_end;
    }

    void init(Signed begin, Signed end, direction dir) {
        if (begin > end) {
            auto buf = begin;
            begin = end + 1;
            end = buf + 1;
        }
        if (dir == direction::forward) {
            m_begin = begin;
            m_end = end;
            m_step = 1;
        } else {
            m_begin = end - 1;
            m_end = begin - 1;
            m_step = -1;
        }
        m_current = m_begin;
    }

    iteration(Signed begin, Signed end, direction dir = direction::forward) {
        init(begin, end, dir);
    }

private:
    Signed m_begin;
    Signed m_end;
    Signed m_current;
    Signed m_step;
};

} // namespace itr
