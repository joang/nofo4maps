# Fixed Point

This project computes the fixed point of a Poincaré map derived from the Hénon-Heiles system, specifically for the case where {x=0} and at a fixed energy level h_0.

## File Description

- **`ode.eq`, `ode.c`, and `ode.h`**: Symbolic links to the model located in the model's folder.
- **`la2d.h`, `la2d.c`**: Linear algebra for 2-by-2 matrices.
- **`newton.c` and `newton.h`**: These files compute a fixed point for the Poincaré map, given a sufficiently accurate initial guess.
- **`main.c`**: Contains the main code that utilizes the provided initial guesses.

## Compiling and Execution

To compile and run the program, follow these steps:

1. Use the provided `Makefile` to compile the code by executing:
   ```bash
   make
   ```
2. Run the executable with an output filename:
   ```bash
   ./main.exe pfix
   ```
