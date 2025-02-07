#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <immintrin.h>

// Simple support class for a 2D vector
class vec2D {
public:
    float x, y;

    // Default constructor initializes both components to 0
    vec2D() { x = y = 0.f; };

    // Constructor initializes components with given values
    vec2D(float _x, float _y) : x(_x), y(_y) {}

    // Constructor initializes components from a vec4
    vec2D(vec4 v) {
        x = v[0];
        y = v[1];
    }

    // Display the vector components
    void display() { std::cout << x << '\t' << y << std::endl; }

    // Overloaded subtraction operator for vector subtraction
    vec2D operator- (vec2D& v) {
        vec2D q;
        q.x = x - v.x;
        q.y = y - v.y;
        return q;
    }
};

// Class representing a triangle for rendering purposes
class triangle {
    Vertex v[3];       // Vertices of the triangle
    float area;        // Area of the triangle
    colour col[3];     // Colors for each vertex of the triangle

public:
    // Constructor initializes the triangle with three vertices
    // Input Variables:
    // - v1, v2, v3: Vertices defining the triangle
    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        // Calculate the 2D area of the triangle
        vec2D e1 = vec2D(v[1].p - v[0].p);
        vec2D e2 = vec2D(v[2].p - v[0].p);
        area = abs(e1.x * e2.y - e1.y * e2.x);
    }

    // Compute edge function
    inline float edgeFunction(const vec2D& a, const vec2D& b, const vec2D& p) const {
        return (p.x - a.x) * (b.y - a.y) - (p.y - a.y) * (b.x - a.x);
    }

    // Helper function to compute the cross product for barycentric coordinates
    // Input Variables:
    // - v1, v2: Edges defining the vector
    // - p: Point for which coordinates are being calculated
    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }

    // Compute barycentric coordinates for a given point
    // Input Variables:
    // - p: Point to check within the triangle
    // Output Variables:
    // - alpha, beta, gamma: Barycentric coordinates of the point
    // Returns true if the point is inside the triangle, false otherwise
    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) / area;
        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) / area;
        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) / area;

        if (alpha < 0.f || beta < 0.f || gamma < 0.f) return false;
        return true;
    }

    // Template function to interpolate values using barycentric coordinates
    // Input Variables:
    // - alpha, beta, gamma: Barycentric coordinates
    // - a1, a2, a3: Values to interpolate
    // Returns the interpolated value
    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }

    // Draw the triangle on the canvas using SIMD
    // Input Variables:
    // - renderer: Renderer object for drawing
    // - L: Light object for shading calculations
    // - ka, kd: Ambient and diffuse lighting coefficients
    void draw(Renderer& renderer, Light& L, float ka, float kd) {
        vec2D minV, maxV;

        // Get the screen-space bounds of the triangle
        getBoundsWindow(renderer.canvas, minV, maxV);

        // Skip very small triangles
        if (area < 1.f) return;

        // Convert vertices to vec2D for easier calculations
        vec2D p0(v[0].p), p1(v[1].p), p2(v[2].p);

        // Compute area in **integer space** for stable barycentric interpolation
        float totalArea = edgeFunction(p0, p1, p2);

        // Precompute edge function values at the top-left corner of the bounding box
        float w0_row = edgeFunction(p1, p2, minV);
        float w1_row = edgeFunction(p2, p0, minV);
        float w2_row = edgeFunction(p0, p1, minV);

        // Compute correct per-pixel increments
        float w0_dx = -(p1.y - p2.y);
        float w0_dy = (p1.x - p2.x);
        float w1_dx = -(p2.y - p0.y);
        float w1_dy = (p2.x - p0.x);
        float w2_dx = -(p0.y - p1.y);
        float w2_dy = (p0.x - p1.x);

        // Iterate over the bounding box and check each pixel
        for (int y = (int)minV.y; y < (int)ceil(maxV.y); y++) {
            float w0 = w0_row;
            float w1 = w1_row;
            float w2 = w2_row;

            for (int x = (int)minV.x; x < (int)ceil(maxV.x); x++) {
                // Check if the pixel is inside the triangle
                if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                    // **Ensure original interpolation order**
                    float alpha = w0 / totalArea;
                    float beta = w1 / totalArea;
                    float gamma = w2 / totalArea;

                    // Interpolate color, depth, and normals exactly as before
                    colour c = interpolate(beta, gamma, alpha, v[0].rgb, v[1].rgb, v[2].rgb);
                    c.clampColour();
                    float depth = interpolate(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
                    vec4 normal = interpolate(beta, gamma, alpha, v[0].normal, v[1].normal, v[2].normal);

                    // **Fix: Ensure normal is normalised before lighting calculation**
                    normal.normalise();

                    // Perform Z-buffer test and apply shading
                    if (renderer.zbuffer(x, y) > depth && depth > 0.01f) {
                        // **Fix: Ensure shading calculation matches the original**
                        L.omega_i.normalise();
                        float dot = max(vec4::dot(L.omega_i, normal), 0.0f);
                        colour a = (c * kd) * (L.L * dot + (L.ambient * kd));

                        unsigned char r, g, b;
                        a.toRGB(r, g, b);
                        renderer.canvas.draw(x, y, r, g, b);
                        renderer.zbuffer(x, y) = depth;
                    }
                }
                // Increment edge function values for next pixel in x direction
                w0 += w0_dx;
                w1 += w1_dx;
                w2 += w2_dx;
            }

            // Increment edge function values for next row in y direction
            w0_row += w0_dy;
            w1_row += w1_dy;
            w2_row += w2_dy;
        }
    }

    // Compute the 2D bounds of the triangle
    // Output Variables:
    // - minV, maxV: Minimum and maximum bounds in 2D space
    void getBounds(vec2D& minV, vec2D& maxV) {
        minV = vec2D(v[0].p);
        maxV = vec2D(v[0].p);
        for (unsigned int i = 1; i < 3; i++) {
            minV.x = min(minV.x, v[i].p[0]);
            minV.y = min(minV.y, v[i].p[1]);
            maxV.x = max(maxV.x, v[i].p[0]);
            maxV.y = max(maxV.y, v[i].p[1]);
        }
    }

    // Compute the 2D bounds of the triangle, clipped to the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    // Output Variables:
    // - minV, maxV: Clipped minimum and maximum bounds
    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        getBounds(minV, maxV);
        minV.x = max(minV.x, 0);
        minV.y = max(minV.y, 0);
        maxV.x = min(maxV.x, canvas.getWidth());
        maxV.y = min(maxV.y, canvas.getHeight());
    }

    // Debugging utility to display the triangle bounds on the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
                canvas.draw(x, y, 255, 0, 0);
            }
        }
    }

    // Debugging utility to display the coordinates of the triangle vertices
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }
};
