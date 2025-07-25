#include <vector>
#include <cmath>
#include <array>
#include "sph.h"

void setVelocity(Particle& p, const std::array<double, 2>& velocity) {
    p.velocity = velocity;
}

void addVelocity(Particle& p, const std::array<double, 2>& acceleration) {
    for (int i = 0; i < 2; ++i) {
        // v = v + a (描画フレームレートでdtは既に考慮済み)
        p.velocity[i] += acceleration[i];
    }
}

void updatePosition(Particle& p) {
    for (int i = 0; i < 2; ++i) { 
        // r = r + v (描画フレームレートでdtは既に考慮済み)
        p.position[i] += p.velocity[i];
    }
}

void init_consts(Particle& p, 
                 double mass,
                 double h,
                 const std::array<double, 2>& position, 
                 const std::array<double, 2>& velocity) {
    p.dt = 0.016; // Default time step (60fps)
    p.mass = mass;
    p.h = h;
    p.position = position;
    p.velocity = velocity;
}

// todo! : 線形探索をグリッド探索に変更する

double ρ (const Particle& p, const std::vector<Particle>& particles) {
    double result = 0.0;
    for (const auto& other : particles) { // neighbor particle search
        if (&other != &p) {
            double dist = std::hypot(other.position[0] - p.position[0], other.position[1] - p.position[1]);
            if (dist < p.h) {
                // ρ = mass * W(r, h)
                result += p.mass * W(dist, p.h);
            }
        }
    }
    return result;
}

// Gaussian Kernel function for SPH
double W(double r, double h) {
    if (r >= 0 && r <= h) {
        double q = r / h;
        // Gaussian kernel: W(r, h) = (1 / (pi * h^2)) * exp(-q^2)
        return (1.0 / (M_PI * h * h)) * exp(-q * q);
    }
    return 0.0;
}

std::array<double, 2> nabla_W(const Particle& p, const Particle& other) {
    std::array<double, 2> result;
    double dist = std::hypot(other.position[0] - p.position[0], other.position[1] - p.position[1]);
    if (dist < p.h) {
        double w = W(dist, p.h);
        // nabla = ∇W(r, h) = (∂W/∂x, ∂W/∂y)
        // nabla = (∂W/∂x, ∂W/∂y) = (∂W/∂r * ∂r/∂x, ∂W/∂r * ∂r/∂y)
        // ∂W/∂r = -2 * W(r, h) / h^2 ・ (x, y)
        // ∂r/∂x = (other.position[0] - p.position[0]) / dist
        // ∂r/∂y = (other.position[1] - p.position[1]) / dist
        result[0] = ((-2 * W(dist, p.h) / (p.h * p.h)) * p.position[0]) * (other.position[0] - p.position[0]) / dist;
        result[1] = ((-2 * W(dist, p.h) / (p.h * p.h)) * p.position[1]) * (other.position[1] - p.position[1]) / dist;
    } else {
        result[0] = 0.0;
        result[1] = 0.0;
    }
    return result;
}

// P is momentum
double P(const Particle& p, const std::vector<Particle>& particles) {
    double result = 0.0;
    for (const auto& other : particles) { // neighbor particle search
        if (&other != &p) {
            double dist = std::hypot(other.position[0] - p.position[0], other.position[1] - p.position[1]);
            if (dist < p.h) {
                // P = mass * W(r, h)
                result += p.mass * W(dist, p.h);
            }
        }
    }
    return result;
}

std::array<double, 2> delta_v(const Particle& p, const std::vector<Particle>& particles) {
    std::array<double, 2> result = {0.0, 0.0};
    for (const auto& other : particles) { // neighbor particle search
        if (&other != &p) {
            double dist = std::hypot(other.position[0] - p.position[0], other.position[1] - p.position[1]);
            if (dist < p.h) {
                // v_n+1 = v_n + (P_i / ρ_i^2) + (P_j / ρ_j^2) *  (∇_i * W(r_i - r_j, h))
                result[0] += (P(p, particles) / (ρ(p, particles) * ρ(p, particles)) + P(other, particles) / (ρ(other, particles) * ρ(other, particles))) * nabla_W(p, other)[0];
                result[1] += (P(p, particles) / (ρ(p, particles) * ρ(p, particles)) + P(other, particles) / (ρ(other, particles) * ρ(other, particles))) * nabla_W(p, other)[1];
            }
        }
    }
    return result;
}

void applyBoundaryConditions(Particle& p, double width, double height) {
    double damping = 0.8; // エネルギー減衰係数
    
    // 左右の境界
    if (p.position[0] < 0) {
        p.position[0] = 0;
        p.velocity[0] = -p.velocity[0] * damping;
    } else if (p.position[0] > width) {
        p.position[0] = width;
        p.velocity[0] = -p.velocity[0] * damping;
    }
    
    // 上下の境界
    if (p.position[1] < 0) {
        p.position[1] = 0;
        p.velocity[1] = -p.velocity[1] * damping;
    } else if (p.position[1] > height) {
        p.position[1] = height;
        p.velocity[1] = -p.velocity[1] * damping;
    }
}

std::array<double, 2> boundaryForce(const Particle& p, double width, double height) {
    std::array<double, 2> force = {0.0, 0.0};
    double boundary_strength = 1.0;
    double boundary_thickness = 10.0;
    
    // 左境界
    if (p.position[0] < boundary_thickness) {
        force[0] += boundary_strength * (boundary_thickness - p.position[0]) / boundary_thickness;
    }
    // 右境界
    if (p.position[0] > width - boundary_thickness) {
        force[0] -= boundary_strength * (p.position[0] - (width - boundary_thickness)) / boundary_thickness;
    }
    // 下境界
    if (p.position[1] < boundary_thickness) {
        force[1] += boundary_strength * (boundary_thickness - p.position[1]) / boundary_thickness;
    }
    // 上境界
    if (p.position[1] > height - boundary_thickness) {
        force[1] -= boundary_strength * (p.position[1] - (height - boundary_thickness)) / boundary_thickness;
    }
    
    return force;
}

std::array<double, 2> centripetal(const Particle& p, double M, double k, int max_x, int max_y) {
    std::array<double, 2> result = {0.0, 0.0};
    std::array<double, 2> center = {static_cast<double>(max_x) / 2.0, static_cast<double>(max_y) / 2.0};

    // dist = ({x, y} - {center_x, center_y})
    // f = -k * M * m / dist^2 (k is a constant, M is the mass of center(virtual point))
    std::array<double, 2> r = {p.position[0] - center[0], p.position[1] - center[1]};
    double dist = std::hypot(r[0], r[1]);
    if (dist > 0) {
        result[0] = -k * M * p.mass * r[0] / (dist * dist);
        result[1] = -k * M * p.mass * r[1] / (dist * dist);
    }

    return result;
}