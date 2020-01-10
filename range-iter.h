// Copyright © 2020 Tokarev Artem. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstdlib>
#include <cmath>


namespace itr {

enum class dir {
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
class range_iter {
public:
    Signed front() const {
        return m_begin;
    }
    Signed back() const {
        return m_end - m_step;
    }
    Signed begin() const {
        return m_begin;
    }
    Signed end() const {
        return m_end;
    }
    Signed next() {
        return m_current + m_step;
    }
    Signed prev() {
        return m_current - m_step;
    }
    bool has_ended() const {
        return m_current == m_end;
    }

    range_iter& operator++() {
        next();
        return *this;
    }
    Signed operator*() {
        return m_current;
    }
    
    void init(Signed begin, Signed end, dir dir) {
        if (begin > end) {
            auto buf = begin;
            begin = end + 1;
            end = buf + 1;
        }
        if (dir == dir::forward) {
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

    range_iter(Signed begin, Signed end, dir dir = dir::forward) {
        init(begin, end, dir);
    }

private:
    Signed m_begin;
    Signed m_end;
    Signed m_current;
    Signed m_step;
};

} // namespace itr
