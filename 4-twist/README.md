# Normal Form

This project computes a twist map from the normal form obtained in **Folder 3**.

## File Descriptions

- **`ode.eq`, `ode.c`, and `ode.h`**: Symbolic links to the model located in the model's folder. These files define the maximum degree for the expansion, with the arithmetic using complex numbers. It suffices to consider half of the symbols and (floor) half of degree that in **Folder 3**.
- **`twist.c` and `twist.h`**: Provide functions for computing the normal form.
- **`util_jet.h`, `util_jet.c`**: Introduces an estimation of validity range of a jet.
- **`main.c`**: Contains the main program logic. It uses the `nofo` function from **Folder 3**.

## Compiling and Execution

To compile and run the program, follow these steps:

1. Use the provided `Makefile` to compile the code by executing:
   ```bash
   make
   ```
2. Run the executable with an output filename:
   ```bash
   ./main.exe nofo twist
   ```
   
## Warnings

The twist can suffer for index selection since only part of the normal form is used, which affects to the twist is obtained. This is user-dependent and it must be studied separately for each problem and interests.
