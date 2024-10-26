# Static Arnold Web Plotting Program

## Software Requirements

- `gfortran` (GNU Fortran compiler)
- `gnuplot` (for plotting)
- `vim` (for editing output files)
- `split` (for splitting files)
- `eps2pdf` or another EPS-to-PDF conversion tool

## Instructions

### Step 1: Generate Initial Data Files

1. Compile and run `arnold.f`:
    ```bash
    gfortran arnold.f
    ./a.out > out
    ```
   - This will produce two files: `out` and `fort.70`.
   - `out`: Contains all data except vertical lines.
   - `fort.70`: Contains data for vertical lines.

### Step 2: Generate Linear Equations

1. Compile and run `linear.f`:
    ```bash
    gfortran linear.f
    ./a.out < out > out2
    ```
   - This will produce two files: `out2` and `fort.71`.
   - `out2`: Contains equations for everything except vertical lines.
   - `fort.71`: Contains the `x` values for vertical lines.
   - Note: The number of lines to be read from `out` and `fort.70` is hard-coded as variables `n` and `m`.

### Step 3: Generate Vertical Lines

1. Compile and run `vertical.f`:
    ```bash
    gfortran vertical.f
    ./a.out
    ```
   - This will produce the file `fort.72`, which contains vertical lines with two extreme `y` values.
   - Note: The number of lines to be read from `fort.71` is hard-coded in the loop over `i`.

### Step 4: Format `out2` for Gnuplot

1. Open `out2` in `vim` and join all lines into a single line:
   - Use `Shift+J` in `vim` to join lines.
2. Save the joined file as `out3`.
3. Add additional details to `out3` to make it a proper `gnuplot` settings file.

### Step 5: Prepare Vertical Line Data

1. Create a folder named `splits` and move `fort.72` there:
    ```bash
    mkdir splits
    cp fort.72 splits/
    ```
2. Inside `splits/`, split `fort.72` into multiple files so each file represents one line:
    ```bash
    cd splits
    split -l <lines_per_file> fort.72
    ```
   - This creates separate files for each vertical line.

3. Create a file named `plot_data` inside `splits/`:
   - Copy all split filenames into `plot_data` and add `w l,` to each line to format for `gnuplot`.
4. Append `plot_data` to `out3` and add other plotting details to finalize the `gnuplot` settings file.

### Step 6: Generate Plot

1. Run `gnuplot` with the settings file:
    ```bash
    gnuplot out3
    ```
   - This will produce an EPS file with omega symbols.

### Step 7: Convert EPS to PDF

1. Convert the EPS file to PDF:
    ```bash
    eps2pdf arnold7.eps arnold7.pdf
    ```

## Notes

Ensure the `n`, `m`, and loop values in `linear.f` and `vertical.f` match your data requirements.
