#include <SFML/Graphics.hpp>
#include <vector>
#include "sph.h"


int main() {
    using screen_size_t = unsigned int;
    screen_size_t width = 800;
    screen_size_t height = 600;
    sf::RenderWindow window(sf::VideoMode(width, height), "Circle List");
    std::vector<sf::CircleShape> circles;

    int gridSize = 7; // Grid size for particle placement
    int totalParticles = gridSize * gridSize; // Total number of particles
    std::vector<Particle> particles(totalParticles);

    // Initialize particles in a grid
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            int index = i * gridSize + j;
            // Initialize each particle with mass, kernel radius, position, and velocity
            init_consts(particles[index],
                        1.0,                                   // mass
                        50.0,                                  // kernel radius
                        std::array<double, 2>{10.0 * i, 10.0 * j},  // position
                        std::array<double, 2>{0.0, 0.0});      // velocity
            updatePosition(particles[index]); // Update position based on velocity and delta time
        }
    }
    

    for (int i = 0; i < totalParticles; ++i) {
        sf::CircleShape c(5.0f); // Circle radius
        c.setFillColor(sf::Color::Cyan);
        c.setPosition(particles[i].position[0], particles[i].position[1]); // Set position based on particle's position
        circles.push_back(c);
    }

    // debug print ρ and P and delta_v and nabla_W
    for (const auto& p : particles) {
        double density = ρ(p, particles);
        double pressure = P(p, particles);
        std::array<double, 2> velocityChange = delta_v(p, particles);
        std::array<double, 2> gradient = nabla_W(p, particles[0]);
        printf("Particle Density: %.3f, Pressure: %.3f, Velocity Change: (%.3f, %.3f), Gradient: (%.3f, %.3f)\n",
               density, pressure, velocityChange[0], velocityChange[1], gradient[0], gradient[1]);
    }

    while (window.isOpen()) {

        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // update velocity
        for (int i = 0; i < totalParticles; ++i) {
            std::array<double, 2> sph_force = delta_v(particles[i], particles);
            std::array<double, 2> boundary_force = boundaryForce(particles[i], width, height);
            
            // 重力も追加
            std::array<double, 2> gravity = {0.0, 9.8 * 0.04}; // 小さくスケール
            
            addVelocity(particles[i], sph_force);
            addVelocity(particles[i], boundary_force);
            addVelocity(particles[i], gravity);
        }
        
        // update positions and apply boundary conditions
        for (int i = 0; i < totalParticles; ++i) {
            updatePosition(particles[i]);
            applyBoundaryConditions(particles[i], width, height); // 境界条件適用
            circles[i].setPosition(particles[i].position[0], particles[i].position[1]);
        }
        
        // update positions of circles based on particles
        for (int i = 0; i < totalParticles; ++i) {
            updatePosition(particles[i]); // Update particle position
            circles[i].setPosition(particles[i].position[0], particles[i].position[1]); // Update circle position
        }

        window.clear();

        for (const auto& c : circles)
            window.draw(c);

        window.display();

        sf::sleep(sf::milliseconds(particles[0].dt * 1000)); // Simulate a frame rate of ~60 FPS
    }

    return 0;
}
