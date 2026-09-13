#pragma once

struct Vec2 {
    double x;
    double y;
};

struct Particle {
    double mass;
    Vec2 position;
    Vec2 velocity;
};
