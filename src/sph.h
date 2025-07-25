#pragma once
#include <array>

struct Particle {
    std::array<double, 2> position; // x, y
    std::array<double, 2> velocity; // vx, vy

    double mass;
    double h; // kernel radius

    double dt = 0.016; // Default time step (60fps)
};

void setVelocity(Particle& p, const std::array<double, 2>& velocity);
void addVelocity(Particle& p, const std::array<double, 2>& acceleration);
void updatePosition(Particle& p);

void init_consts(Particle& p, 
                 double mass,
                 double h,
                 const std::array<double, 2>& position, 
                 const std::array<double, 2>& velocity);

// todo! : 線形探索をグリッド探索に変更する
double W(double r, double h);
std::array<double, 2> nabla_W(const Particle& p, const Particle& other);
double ρ(const Particle& p, const std::vector<Particle>& particles);
double P(const Particle& p, const std::vector<Particle>& particles);
std::array<double, 2> delta_v(const Particle& p, const std::vector<Particle>& particles);
std::array<double, 2> centripetal(const Particle& p, double M, double k, int max_x, int max_y);
void applyBoundaryConditions(Particle& p, double width, double height);
std::array<double, 2> boundaryForce(const Particle& p, double width, double height);
std::array<double, 2> centripetal(const Particle& p, double M, double k, int max_x, int max_y);