#include "io.hpp"

void write_state(
    std::ostream& output,
    double time,
    const std::vector<Particle>& particles){
    for(std::size_t id = 0; id < particles.size(); ++id){
        output << time << ','
               << id << ','
               << particles[id].position.x << ','
               << particles[id].position.y << '\n';
    }
}
