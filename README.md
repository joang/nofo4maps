# Explicit Numerical Computation of Normal Forms for Poincaré Maps

## Introduction

This project provides a basic code sample for computing normal forms and when the the fixed points is elliptic Hamiltonian also twist maps. The sample is provided for a specific model. While the normal form construction does not necessarily require Hamiltonian dynamics, we illustrate the procedure using the classical **Hénon-Heiles** system, a standard Hamiltonian model in celestial mechanics.

For further details, see the manuscript [GJJZ24](http://www.maia.ub.es/dsg/2024/2405gimeno.pdf).

## Requirements

The code is written in C and has been successfully compiled with the GNU C Compiler version 12.2.0. The programs are self-contained, requiring no extra libraries, and they depend on appropriate initial guesses to function correctly.


Output files from the [Taylor v2.2](https://github.com/joang/taylor2-dist) package are provided in the **`0-model`** folder. These include structures like `MY_JET` and `MY_FLOAT`, and various functions. The folder also contains instructions on generating these output files using `taylor` if users wish to extend the sample.

## Usage

To use the code, follow the folder structure in the specified order, adhering to the instructions in the accompanying `README` files:

- **`0-model`**: Contains the model and outputs from `taylor` used in subsequent folders. Everything required for the sample is included here.
- **`1-fixed-point`**: Computes a fixed point of the Poincaré map for the Hénon-Heiles system, fixing an energy level and a spatial section at {x = 0}.
- **`2-jet-fixed-point`**: Computes the high-order derivatives of the Poincaré map at the fixed point computed in `1`.
- **`3-normal-form`**: Performs a change of coordinates on the jet obtained in `2` to derive its normal form.
- **`4-twist`**: Computes a twist map based on the output of `3`, leveraging the elliptic fixed point from `1` and the Hamiltonian nature of the system.


## Warnings

This code is a sample designed to illustrate the process of computing normal forms for Poincaré maps. Specifically, it addresses 2-by-2 systems, solving linear systems and computing eigenpairs.

Note that this sample does not include components such as parameter-dependent normal forms or torus visualizations, as these are outside the scope of the demonstration and might detract from clarity.

## Enhancements

You may consider the following improvements:

- Implement parameter-dependence normal form.
- Develop the torus visualization for a radius within the validity range.
- Improve the linear algebra parts by incorporating custom libraries.
- Run the code with extended precision (this requires modifying the output files from `taylor` in the **`0-model`** folder).
- Use wrappers to run the code, particularly the part of `taylor`, in a different programming language, e.g., Python.

## Authors and Acknowledgments

The authors of the paper and this code project are:

- J. Gimeno
- À. Jorba
- M. Jorba-Cuscó
- M. Zou

**Date**: December 2024.

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

For more information about the GPL v3.0, please visit the official GNU website: [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).
