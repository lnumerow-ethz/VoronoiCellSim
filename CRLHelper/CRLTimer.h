#pragma once

#include <chrono>
#include <iostream>

class CRLTimer {
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

   public:
    void printTime(std::string note) {
        auto t = std::chrono::high_resolution_clock::now();
        int total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t - t0).count();
        int this_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t - t1).count();
        t1 = t;

        std::cout << note << ": " << this_ms << " ms, Total: " << total_ms << " ms" << std::endl;
    }
};
