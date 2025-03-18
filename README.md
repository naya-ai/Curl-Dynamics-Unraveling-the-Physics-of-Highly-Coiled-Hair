# Thesis Research: Curl Dynamics – Unraveling the Physics of Highly Coiled Hair

This repository contains my thesis research project for the Thomas F. Freeman Honors College at Texas Southern University. The work focuses on the realistic simulation of highly coiled hair (Type 4 hair) using a physics-based model.

## Research Motivation

Highly coiled hair has often been misrepresented or oversimplified in computer graphics, which can lead to cultural misunderstandings and a lack of authentic representation. The goal of this research is to:
- Develop a simulation that accurately captures the complex structure and dynamic behavior of highly coiled hair.
- Address the technical challenges of simulating natural hair dynamics, such as stretching, bending, twisting, and the stabilization of twist via phase locking.
- Contribute to more inclusive digital media by providing a model that better represents the diversity of hair textures, specifically Type 4 hair.

## What This Project Does

This project uses a 3D Discrete Elastic Rod (DER) approach to simulate hair strands. Key aspects include:

- **Node-Based Representation:**  
  The hair strand is discretized into multiple connected nodes, with each segment having its own rest length, twist angle, and local frame (tangent, normal, and binormal).

- **Physical Forces:**  
  The simulation incorporates:
  - **Stretching:** To maintain natural segment lengths.
  - **Bending:** To simulate curvature along the strand.
  - **Twisting:** To emulate natural hair rotation, with twist angles being wrapped to prevent uncontrolled accumulation.
  - **Phase Locking:** A mechanism that aligns adjacent twist angles, stabilizing the overall curl pattern and reducing erratic rotations.

- **Interactive Parameter Control:**  
  An integrated ImGui-based GUI allows real-time adjustments of:
  - Number of nodes
  - Stretch stiffness (StretchK)
  - Bend stiffness (BendK)
  - Twist stiffness (TwistK)
  - Phase lock stiffness (PhaseLockK)
  - Damping  
  This interactive approach provides immediate visual feedback, enabling a thorough exploration of how each parameter affects the behavior and appearance of highly coiled hair.

## Purpose and Impact

The primary aim of this research is to:
- Demonstrate the feasibility of accurately modeling the complex dynamics of highly coiled hair.
- Provide a tool that can be used to analyze and visualize the influence of different physical parameters on hair behavior.
- Enhance digital representations of Type 4 hair in applications such as video games, animation, and film, thus promoting inclusivity and cultural authenticity in visual media.
