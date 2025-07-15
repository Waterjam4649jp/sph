#pragma once
#include <array>

struct Particle {
    std::array<float, 2> position; // x, y
    std::array<float, 2> velocity; // vx, vy
    float mass;
    float dt;
};

void set_delta_time(Particle& p, float dt);
void setVelocity(Particle& p, const std::array<float, 2>& velocity);
void updatePosition(Particle& p);
