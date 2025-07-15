#include <array>
#include "sph.h"

void set_delta_time(Particle& p, float dt) {
    p.dt = dt;
}

void setVelocity(Particle& p, const std::array<float, 2>& velocity) {
    p.velocity = velocity;
}

void updatePosition(Particle& p) {
    for (int i = 0; i < 2; ++i) {
        p.position[i] += p.velocity[i] * p.dt;
    }
}