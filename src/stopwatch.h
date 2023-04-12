#pragma once

#include <chrono>

class Stopwatch
{
public:
	void start() { m_start = std::chrono::system_clock::now(); }
	void stop() { m_end = std::chrono::system_clock::now(); }
	float getMs() { return std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_start).count() / 1000.f / 1000.f; }

private:
	std::chrono::time_point<std::chrono::system_clock> m_start, m_end;
};
