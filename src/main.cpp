#include <SFML/Graphics.hpp>
#include <vector>
#include "sph.h"

int main() {
    sf::RenderWindow window(sf::VideoMode(800, 600), "Circle List");

    std::vector<sf::CircleShape> circles;

    int numParticles = 5; // Number of circles to create
    Particle particles[numParticles];

    for (int i = 0; i < numParticles; ++i) {
        particles[i].position = {static_cast<float>(i * 60 + 50), 100.0f}; // x, y
        particles[i].velocity = {20.0f, 20.0f}; // vx, vy
        particles[i].mass = 1.0f; // mass

        set_delta_time(particles[i], 0.016f); // Set delta time to ~60 FPS
        updatePosition(particles[i]); // Update position based on velocity and delta time
    }

    for (int i = 0; i < numParticles; ++i) {
        sf::CircleShape c(30);
        c.setFillColor(sf::Color::Cyan);
        c.setPosition(particles[i].position[0], particles[i].position[1]); // Set position based on particle's position
        circles.push_back(c);
    }

    while (window.isOpen()) {

        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        
        // udpate positions of circles based on particles
        for (int i = 0; i < numParticles; ++i) {
            updatePosition(particles[i]); // Update particle position
            circles[i].setPosition(particles[i].position[0], particles[i].position[1]); // Update circle position
        }

        window.clear();

        for (const auto& c : circles)
            window.draw(c);

        window.display();
    }

    return 0;
}
