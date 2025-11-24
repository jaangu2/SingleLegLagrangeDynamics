# SingleLegLagrangeDynamics
This project runs numerical simulation of a two-segment robot leg dynamically moving through stance phase.

## Overview:
The pdf includes the hand derivations of the equations defining the two-segment leg system with a spring on the hip and knee joints. \
The FullStanceSolverWithVaryingK1.m file runs the numerical simulation and plots the results. \
The solveForK1.m file is a helper file which solves for the all-important K1 value for each timestep of the simulation.

## Purpose:
The point of this project is to prove the concept that adjusting the hip stiffness through changing the K1 value allows for ground reaction force redirection. Once this concept is proven, moving onto a full quadrupedal system is the next step.

## Method:
I used symbolic variables to set up my dynamic equations using the Lagrange method. By setting up the potential and kinetic energy equations, and including nonconservative forces, the Lagrange equation can be used to find the dynamic equations of this system, which are then solved using ode45 for every timestep through a single dynamic stance phase.
