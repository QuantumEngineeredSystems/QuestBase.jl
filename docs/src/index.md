# The Common Interface for the QUEST ecosystem

The QUEST common interface ties together the utilities of the ecosystem into a single unified interface.
It is designed ensure compatiabilty between the different available packages.

This documentation is made to pool together the docs of the various Quest libraries
to paint the overarching picture and document the shared/common functionality.

## Packages of Quest Ecosystem

### [HarmonicBalance.jl](https://github.com/QuantumEngineeredSystems/HarmonicBalance.jl):

A package for applying the harmonic balance method to classical driven nonlinear resonators.
It computes the stroboscopic effective equations of motion of the system at the characteristic response frequencies of the system.
Both [Krylov-Bogoliubov](https://quantumengineeredsystems.github.io/HarmonicBalance.jl/stable/manual/extracting_harmonics#Krylov-Bogoliubov) averaging method to higher orders and the [harmonic balance method](https://quantumengineeredsystems.github.io/HarmonicBalance.jl/stable/manual/extracting_harmonics#Krylov-Bogoliubov) are implemented.

### [HarmonicSteadyState.jl](https://github.com/QuantumEngineeredSystems/HarmonicSteadyState.jl)

A package for computing the classical steady state of the effective stroboscopic systems driven nonlinear resonators. Given one has the autonomous equations of motion of the system in the rotating frame of the characteristic response frequencies, it collect steady states methods to find and describe the stationary responses of the system. It supports the following methods:

- fixed point steady states with Homotopy Continuation
- Finding Limit-cycle  with Homotopy Continuation
- Stability analysis
- Linear response of the steady state in the (non-)rotating frame
- Parameter sweeps
- Plotting utilities
