/*
    DiscreteElasticRod3D.cpp

    Overview:
      This file implements a physics-based simulation to model highly coiled hair (Type 4 hair) 
      using a 3D Discrete Elastic Rod (DER) approach. The simulation initializes a tight-coiled rod 
      that mimics the complex geometry of highly coiled hair and applies various physical forces to 
      capture realistic dynamics.

    Key Features:
      • Node-Based Representation:
          The rod is discretized into a series of nodes, with each segment having associated rest 
          lengths, twist angles, and local frames (tangent, normal, and binormal).

      • Physical Forces Modeled:
          - Stretching: Ensures each rod segment maintains its natural length.
          - Bending: Simulates curvature by applying torques based on the angle between adjacent 
            segments.
          - Twisting: Introduces twist along the rod’s axis; twist angles are wrapped to prevent 
            uncontrolled rotations.
          - Phase Locking: A novel mechanism that penalizes differences between adjacent twist 
            angles, helping to maintain a stable and consistent curl pattern.

      • Advanced Spectral Smoothing (Inspired by Curly-Cue):
          A new smoothing method is applied to the twist angles to mimic a Fourier-based phase locking
          mechanism. (Currently implemented as a moving-average filter, with the option to replace it 
          with an FFT-based method for more precise spectral analysis.)

      • Interactive Parameter Control:
          An ImGui-based GUI provides real-time sliders for:
            - Number of nodes
            - Stretch stiffness (StretchK)
            - Bend stiffness (BendK)
            - Twist stiffness (TwistK)
            - Phase lock stiffness (PhaseLockK)
            - Damping factor
            - Spectral Smoothing toggle and window size
          This enables immediate visual feedback on how each parameter affects the simulation.

      • Rendering and Interaction:
          Utilizes GLFW for windowing and input, and OpenGL for rendering the simulation. The code 
          includes an adjustable camera system, allowing for dynamic visualization of the rod's behavior.

    Purpose:
      The simulation serves as a proof-of-concept to demonstrate the feasibility of accurately modeling 
      the dynamics of highly coiled hair. By adjusting the parameters, the model can represent a range 
      of behaviors—from natural, stable curls to dynamic, bouncy movements. This work not only addresses 
      technical challenges in hair simulation but also contributes to the inclusive representation of Type 4 
      hair in digital media, with potential applications in video games, animation, and film.

    Usage:
      - Compile and run the project to visualize the simulation.
      - Use the provided GUI sliders to experiment with different parameter combinations and observe their 
        effects on the hair simulation.
      - A screen recording demonstrating various parameter combinations will be attached for further review.

    Dependencies:
      - GLFW for window management and input handling.
      - OpenGL for rendering graphics.
      - ImGui for GUI creation and interactive control.
      - Standard C++ libraries for mathematical operations and data structures.
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

#ifdef _WIN32
#include <windows.h>
#endif

#include <GLFW/glfw3.h>
#define GL_SILENCE_DEPRECATION
#include <GL/glu.h>

// ImGui includes
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// --------------------------------------------------
// Basic Vec3 & Math Helpers
// --------------------------------------------------
struct Vec3 {
    float x, y, z;
};

inline Vec3 operator+(const Vec3 &a, const Vec3 &b) { 
    return {a.x + b.x, a.y + b.y, a.z + b.z}; 
}

inline Vec3 operator-(const Vec3 &a, const Vec3 &b) { 
    return {a.x - b.x, a.y - b.y, a.z - b.z}; 
}

inline Vec3 operator*(const Vec3 &v, float s) { 
    return {v.x * s, v.y * s, v.z * s}; 
}

inline Vec3 operator*(float s, const Vec3 &v) { 
    return {v.x * s, v.y * s, v.z * s}; 
}

inline Vec3 operator-(const Vec3 &v) { 
    return {-v.x, -v.y, -v.z}; 
}

// Compound assignment operators for Vec3
inline Vec3& operator+=(Vec3 &a, const Vec3 &b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

inline Vec3& operator-=(Vec3 &a, const Vec3 &b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

inline float dot(const Vec3 &a, const Vec3 &b) { 
    return a.x * b.x + a.y * b.y + a.z * b.z; 
}

inline float length(const Vec3 &v) { 
    return std::sqrt(dot(v, v)); 
}

inline Vec3 normalize(const Vec3 &v) {
    float l = length(v);
    if (l < 1e-8f) return {0.f, 0.f, 0.f};
    return {v.x / l, v.y / l, v.z / l};
}

inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
}

// Custom clamp function
inline float clampf(float val, float minVal, float maxVal) {
    if (val < minVal) return minVal;
    if (val > maxVal) return maxVal;
    return val;
}

// Wrap angle to [-PI, PI]
inline float wrapAngle(float angle) {
    const float PI = 3.14159f;
    const float TWO_PI = 6.28318f;
    while (angle > PI) angle -= TWO_PI;
    while (angle < -PI) angle += TWO_PI;
    return angle;
}

// --------------------------------------------------
// Frame3: holds tangent, normal, binormal
// --------------------------------------------------
struct Frame3 {
    Vec3 tangent;
    Vec3 normal;
    Vec3 binormal;
};

// --------------------------------------------------
// Discrete Elastic Rod Class
// --------------------------------------------------
class DiscreteElasticRod {
public:
    std::vector<Vec3> positions;
    std::vector<Vec3> velocities;
    std::vector<Frame3> frames;
    std::vector<float> restLengths;
    std::vector<float> restTwists;
    std::vector<float> twistAngles;

    // Physical parameters
    float stretchK  = 100.f;
    float bendK     = 5.f;
    float twistK    = 0.5f;
    float phaseLockK= 0.f;
    float mass      = 0.01f;
    float damping   = 0.995f;
    float timeStep  = 0.005f;
    bool  fixFirst  = true;
    Vec3  gravity   = {0.f, -9.8f, 0.f};

    // New parameters for spectral smoothing (placeholder for FFT-based filtering)
    bool useSpectralSmoothing = false;
    int smoothingWindow = 3; // should be odd for symmetric smoothing

    DiscreteElasticRod() {}

    // Initialize a coiled rod with a tight-coil configuration
    void initRod(int numNodes) {
        positions.resize(numNodes);
        velocities.resize(numNodes, {0.f, 0.f, 0.f});
        frames.resize(numNodes - 1);
        restLengths.resize(numNodes - 1);
        restTwists.resize(numNodes - 1, 0.f);
        twistAngles.resize(numNodes - 1, 0.f);

        float coilRadius = 0.05f;
        float coilPitch  = 0.02f;
        float coilTurns  = 5.0f;
        int nSeg = numNodes - 1;

        for (int i = 0; i < numNodes; i++) {
            float t = (float)i / (float)nSeg * coilTurns * 2.f * 3.14159f;
            float y = coilPitch * i;
            float x = coilRadius * std::cos(t);
            float z = coilRadius * std::sin(t);
            positions[i] = {x, y, z};
        }
        for (int i = 0; i < nSeg; i++) {
            Vec3 d = positions[i + 1] - positions[i];
            restLengths[i] = length(d);
        }
        updateFrames();
    }

    // Update the local frame for each segment
    void updateFrames() {
        int nSeg = (int)frames.size();
        for (int i = 0; i < nSeg; i++) {
            Vec3 edge = positions[i + 1] - positions[i];
            float L = length(edge);
            Frame3 &f = frames[i];
            if (L < 1e-8f) {
                f.tangent  = {0.f, 1.f, 0.f};
                f.normal   = {1.f, 0.f, 0.f};
                f.binormal = {0.f, 0.f, 1.f};
                continue;
            }
            f.tangent = normalize(edge);
            if (i == 0) {
                Vec3 guess = {0.f, 1.f, 0.f};
                if (std::fabs(dot(guess, f.tangent)) > 0.99f)
                    guess = {1.f, 0.f, 0.f};
                Vec3 n = normalize(cross(f.tangent, guess));
                Vec3 b = cross(f.tangent, n);
                f.normal = n;
                f.binormal = b;
            } else {
                Frame3 &prev = frames[i - 1];
                Vec3 projN = prev.normal - f.tangent * dot(prev.normal, f.tangent);
                float lenN = length(projN);
                if (lenN < 1e-8f) {
                    f.normal = {1.f, 0.f, 0.f};
                    f.binormal = cross(f.tangent, f.normal);
                } else {
                    projN = projN * (1.f / lenN);
                    Vec3 projB = prev.binormal - f.tangent * dot(prev.binormal, f.tangent);
                    projB = normalize(projB);
                    float alpha = twistAngles[i];
                    f.normal = rotateAroundAxis(projN, f.tangent, alpha);
                    f.binormal = rotateAroundAxis(projB, f.tangent, alpha);
                }
            }
        }
    }

    // Rotate vector v around axis by angle alpha
    Vec3 rotateAroundAxis(const Vec3 &v, const Vec3 &axis, float alpha) {
        Vec3 a = normalize(axis);
        float c = std::cos(alpha);
        float s = std::sin(alpha);
        Vec3 axv = cross(a, v);
        float adv = dot(a, v);
        return v * c + axv * s + a * (adv * (1.f - c));
    }

    // Compute forces and twist rates acting on the rod
    void computeForces(std::vector<Vec3> &outF, std::vector<float> &outTwistRates) {
        int n = (int)positions.size();
        int nSeg = n - 1;
        // Gravity
        for (int i = 0; i < n; i++)
            outF[i] = {mass * gravity.x, mass * gravity.y, mass * gravity.z};
        // Reset twist rates
        for (int i = 0; i < nSeg; i++)
            outTwistRates[i] = 0.f;
        // Stretching
        for (int i = 0; i < nSeg; i++) {
            Vec3 d = positions[i + 1] - positions[i];
            float L = length(d);
            if (L < 1e-8f) continue;
            float stretch = L - restLengths[i];
            float Fmag = stretchK * stretch;
            Vec3 dir = (1.f / L) * d;
            outF[i] += dir * Fmag;
            outF[i + 1] -= dir * Fmag;
        }
        // Bending
        for (int i = 1; i < nSeg; i++) {
            Vec3 t1 = frames[i - 1].tangent;
            Vec3 t2 = frames[i].tangent;
            float dotT = clampf(dot(t1, t2), -1.f, 1.f);
            float angle = std::acos(dotT);
            float bendEnergy = bendK * angle;
            Vec3 axis = cross(t1, t2);
            float lenA = length(axis);
            if (lenA < 1e-8f) continue;
            Vec3 axisDir = (1.f / lenA) * axis;
            float torqueMag = bendEnergy;
            Vec3 force_i = cross(axisDir, t1) * torqueMag;
            Vec3 force_i1 = -cross(axisDir, t2) * torqueMag;
            outF[i] += force_i;
            outF[i + 1] += force_i1;
        }
        // Twisting
        for (int i = 0; i < nSeg; i++) {
            float twistDiff = twistAngles[i] - restTwists[i];
            float torque = twistK * twistDiff;
            outTwistRates[i] = -torque;
        }
        // Phase locking to stabilize twist differences
        if (phaseLockK > 0.f) {
            for (int i = 1; i < nSeg; i++) {
                float phaseDiff = twistAngles[i] - twistAngles[i - 1];
                float lockTorque = phaseLockK * phaseDiff;
                outTwistRates[i] -= lockTorque;
                outTwistRates[i - 1] += lockTorque;
            }
        }
    }

    // Apply a simple moving-average filter to twist angles
    void smoothTwistAngles() {
        int nSeg = twistAngles.size();
        if (nSeg < smoothingWindow) return;
        std::vector<float> smoothed(nSeg, 0.f);
        int halfWindow = smoothingWindow / 2;
        for (int i = 0; i < nSeg; i++) {
            float sum = 0.f;
            int count = 0;
            for (int j = i - halfWindow; j <= i + halfWindow; j++) {
                if (j >= 0 && j < nSeg) {
                    sum += twistAngles[j];
                    count++;
                }
            }
            smoothed[i] = sum / count;
        }
        twistAngles = smoothed;
    }

    // Integrate positions and twist angles using Euler integration
    void integrate() {
        int n = (int)positions.size();
        int nSeg = n - 1;
        std::vector<Vec3> forces(n, {0.f, 0.f, 0.f});
        std::vector<float> twistRates(nSeg, 0.f);
        computeForces(forces, twistRates);
        // Integrate positions
        for (int i = 0; i < n; i++) {
            if (fixFirst && i == 0) {
                velocities[i] = {0.f, 0.f, 0.f};
                continue;
            }
            Vec3 accel = {forces[i].x / mass, forces[i].y / mass, forces[i].z / mass};
            velocities[i] = velocities[i] + accel * timeStep;
            positions[i] = positions[i] + velocities[i] * timeStep;
        }
        // Integrate twist angles and wrap them
        for (int i = 0; i < nSeg; i++) {
            float angVel = twistRates[i] / mass;
            twistAngles[i] += angVel * timeStep;
            twistAngles[i] = wrapAngle(twistAngles[i]);
        }
        // Optionally apply spectral (moving-average) smoothing
        if (useSpectralSmoothing)
            smoothTwistAngles();
        // Apply damping to velocities
        for (auto &v : velocities)
            v = v * damping;
        updateFrames();
    }

    // Perform one simulation step
    void step() {
        integrate();
    }
};

// --------------------------------------------------
// GLOBALS
// --------------------------------------------------
static DiscreteElasticRod gRod;
static int   gNumNodes = 20;
static float gCamDist  = 2.5f;
static float gCamYaw   = 0.f;
static float gCamPitch = 0.5f;
static bool  gLeftDown = false;
static double gPrevMouseX = 0.0;
static double gPrevMouseY = 0.0;

// Mouse and camera callbacks
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if(button == GLFW_MOUSE_BUTTON_LEFT) {
        if(action == GLFW_PRESS) {
            gLeftDown = true;
            glfwGetCursorPos(window, &gPrevMouseX, &gPrevMouseY);
        } else if(action == GLFW_RELEASE) {
            gLeftDown = false;
        }
    }
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    gCamDist *= std::exp(-0.1f * (float)yoffset);
    if(gCamDist < 0.1f)  gCamDist = 0.1f;
    if(gCamDist > 20.f)  gCamDist = 20.f;
}

void updateCamera(GLFWwindow* window) {
    if(gLeftDown) {
        double mx, my;
        glfwGetCursorPos(window, &mx, &my);
        float dx = (float)(mx - gPrevMouseX) * 0.005f;
        float dy = (float)(my - gPrevMouseY) * 0.005f;
        gPrevMouseX = mx;
        gPrevMouseY = my;
        gCamYaw   += dx;
        gCamPitch += dy;
        if(gCamPitch > 1.5f)  gCamPitch = 1.5f;
        if(gCamPitch < -1.5f) gCamPitch = -1.5f;
    }
}

void setCamera(int w, int h) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = (float)w / (float)h;
    gluPerspective(45.f, aspect, 0.01f, 100.f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    float cx = gCamDist * std::cos(gCamPitch) * std::sin(gCamYaw);
    float cy = gCamDist * std::sin(gCamPitch);
    float cz = gCamDist * std::cos(gCamPitch) * std::cos(gCamYaw);
    gluLookAt(cx, cy, cz,  0.f, 0.f, 0.f,  0.f, 1.f, 0.f);
}

// --------------------------------------------------
// Render the rod and its local frames
// --------------------------------------------------
void renderRod(const DiscreteElasticRod &rod) {
    int n = (int)rod.positions.size();
    if(n < 2) return;
    glLineWidth(2.f);
    glColor3f(1.f, 1.f, 1.f);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < n; i++) {
        const Vec3 &p = rod.positions[i];
        glVertex3f(p.x, p.y, p.z);
    }
    glEnd();
    glPointSize(6.f);
    glBegin(GL_POINTS);
    for (int i = 0; i < n; i++) {
        if (i == 0 && rod.fixFirst)
            glColor3f(1.f, 1.f, 0.f);
        else
            glColor3f(1.f, 1.f, 1.f);
        const Vec3 &p = rod.positions[i];
        glVertex3f(p.x, p.y, p.z);
    }
    glEnd();
    // Render local frames
    for (int i = 0; i < n - 1; i++) {
        const Frame3 &f = rod.frames[i];
        Vec3 mid = (rod.positions[i] + rod.positions[i + 1]) * 0.5f;
        float scale = 0.05f;
        // Tangent (red)
        Vec3 tEnd = mid + f.tangent * scale;
        glColor3f(1.f, 0.f, 0.f);
        glBegin(GL_LINES);
        glVertex3f(mid.x, mid.y, mid.z);
        glVertex3f(tEnd.x, tEnd.y, tEnd.z);
        glEnd();
        // Normal (green)
        Vec3 nEnd = mid + f.normal * scale;
        glColor3f(0.f, 1.f, 0.f);
        glBegin(GL_LINES);
        glVertex3f(mid.x, mid.y, mid.z);
        glVertex3f(nEnd.x, nEnd.y, nEnd.z);
        glEnd();
        // Binormal (blue)
        Vec3 bEnd = mid + f.binormal * scale;
        glColor3f(0.f, 0.f, 1.f);
        glBegin(GL_LINES);
        glVertex3f(mid.x, mid.y, mid.z);
        glVertex3f(bEnd.x, bEnd.y, bEnd.z);
        glEnd();
    }
}

// --------------------------------------------------
// MAIN
// --------------------------------------------------
int main() {
    if(!glfwInit()){
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }
    GLFWwindow* window = glfwCreateWindow(1280, 720, "3D DER - Enhanced Phase Locking", NULL, NULL);
    if(!window){
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetScrollCallback(window, scrollCallback);

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    const char* glsl_version = "#version 130";
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Initialize rod
    gRod.initRod(gNumNodes);

    // Main loop
    while(!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        updateCamera(window);
        gRod.step();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.15f, 0.15f, 0.2f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        setCamera(w, h);
        renderRod(gRod);

        // ImGui UI
        ImGui::Begin("DER Hair Controls");
        if(ImGui::SliderInt("Num Nodes", &gNumNodes, 5, 50))
            gRod.initRod(gNumNodes);
        ImGui::SliderFloat("StretchK", &gRod.stretchK, 0.f, 300.f);
        ImGui::SliderFloat("BendK", &gRod.bendK, 0.f, 50.f);
        ImGui::SliderFloat("TwistK", &gRod.twistK, 0.f, 5.f);
        ImGui::SliderFloat("PhaseLockK", &gRod.phaseLockK, 0.f, 10.f, "Lock=%.3f");
        ImGui::SliderFloat("Damping", &gRod.damping, 0.9f, 1.f);
        ImGui::Checkbox("Fix First Node", &gRod.fixFirst);
        ImGui::Checkbox("Spectral Smoothing", &gRod.useSpectralSmoothing);
        ImGui::SliderInt("Smoothing Window", &gRod.smoothingWindow, 3, 15);
        if(ImGui::Button("Reset Rod"))
            gRod.initRod(gNumNodes);
        ImGui::Text("Camera Dist=%.2f Yaw=%.2f Pitch=%.2f", gCamDist, gCamYaw, gCamPitch);
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
