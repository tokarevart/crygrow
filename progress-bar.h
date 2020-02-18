// Copyright © 2020 Artyom Tokarev. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <iostream>
#include <string>

class progress_bar {
public:
    void restart() {
        m_count = 0;
        display();
    }
    void stop() {
        display();
        std::cout << std::endl;
        m_is_started = false;
    }

    std::string header() const {
        return m_header;
    }
    std::size_t bar_width() const {
        return m_bar_width;
    }
    std::size_t max() const {
        return m_max;
    }
    std::size_t count() const {
        return m_count;
    }
    // return new count value
    std::size_t set_count(std::size_t new_value) {
        if (new_value < max()) {
            m_count = new_value;
            display();
        } else {
            m_count = max();
            stop();
        }
        return m_count;
    }

    std::size_t operator+=(std::size_t increment) {
        return set_count(count() + increment);
    }
    std::size_t operator-=(std::size_t decrement) {
        return set_count(count() - decrement);
    }
    std::size_t operator++() {
        return set_count(count() + 1);
    }
    std::size_t operator--() {
        return set_count(count() - 1);
    }

    progress_bar(std::string header, std::size_t max, std::size_t width,
                 char complete = '=', char incomplete = ' ')
        : m_header(std::move(header)), m_max(max), m_bar_width(width),
        m_complete_char(complete), m_incomplete_char(incomplete) {}


private:
    bool m_is_started = true;
    std::size_t m_count = 0;

    const std::string m_header;
    const std::size_t m_max;
    const std::size_t m_bar_width;
    const char m_complete_char;
    const char m_incomplete_char;

    void display() const {
        if (!m_is_started)
            return;

        auto completed_width = (bar_width() * count()) / max();

        std::string bar = "\r" + header() + " [";
        bar.reserve(header().size() + bar_width() + 12);
        for (std::size_t i = 0; i < completed_width; ++i)
            bar += m_complete_char;
        for (std::size_t i = completed_width; i < bar_width(); ++i)
            bar += m_incomplete_char;
        bar += "] " + std::to_string((count() * 100) / max()) + "%";

        std::cout << bar;
        std::cout.flush();
    }
};
