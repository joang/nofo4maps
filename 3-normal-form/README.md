# Normal Form

This project computes the normal form the local approximation obtained in **Folder 2**.

## File Descriptions

- **`ode.eq`, `ode.c`, and `ode.h`**: Symbolic links to the model located in the model's folder. These files define the maximum degree for the expansion, with the arithmetic using complex numbers.
- **`la2d.h`, `la2d.c`**: Implementations for linear algebra operations on 2-by-2 matrices.
- **`nofo.c` and `nofo.h`**: Provide functions for computing the normal form.
- **`util_jet.h`, `util_jet.c`**: Introduces an estimation of validity range of a jet.
- **`main.c`**: Contains the main program logic. It uses the `jet_pfix` function from **Folder 2**.

## Compiling and Execution

To compile and run the program, follow these steps:

1. Use the provided `Makefile` to compile the code by executing:
   ```bash
   make
   ```
2. Run the executable with an output filename:
   ```bash
   ./main.exe jet_pfix nofo
   ```
   
## Warnings

Some parts are model dependent such as the function `is_unavoidable_resonant` in `nofo` files. 
The function `is_resonant_monomial` only detects possible unavoidable resonances for a given log10 tolerance and it may need to be adapted for different scenarios.
