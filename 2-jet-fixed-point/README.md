# Jet at the Fixed Point

This project computes high-order derivatives the Poincaré map at the fixed point of the Hénon-Heiles system with a fixed energy level h_0 and spatial section {x=0}.

## File Description

- **`ode.eq`, `ode.c`, and `ode.h`**: Symbolic links to the model located in the model's folder. It sets the maximum degree for the expansion
- **`poinca.c` and `poinca.h`**: These files construct the high-order derivatives of the Poincaré map and are model-dependent.
- **`main.c`**: Contains the main code. It uses the `pfix` from the 1 folder.

## Compiling and Execution

To compile and run the program, follow these steps:

1. Use the provided `Makefile` to compile the code by executing:
   ```bash
   make
   ```
2. Run the executable with an output filename:
   ```bash
   ./main.exe pfix jet_pfix
   ```
