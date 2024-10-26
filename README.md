# Wavelet Transform for Analyzing Molecular Simulations

This repository contains a Fortran program for wavelet transform analysis, which enables the study of Intramolecular Vibrational Redistribution (IVR) in high-dimensional molecular simulation trajectories. By using dimensionality reduction techniques with the wavelet transform, we can understand the complex problem of IVR in molecular systems in terms of the dominant atomic motions, making the analysis more tractable.

## Background
The wavelet transform helps analyze vibrational frequencies and time-dependent dynamics in molecular simulations, assisting in the study of energy transfer between the dominant vibrational modes of the molecule.

For a thorough understanding of wavelet transform methods applied here, refer to:
1. Vela-Arevalo LV, Wiggins S. *Time-frequency analysis of classical trajectories of polyatomic molecules.* International Journal of Bifurcation and Chaos. 2001; 11(05):1359-80.
2. Rahaman A, Wheeler RA. *Wavelet transforms for determining time-dependent vibrational frequencies.* Journal of Chemical Theory and Computation. 2005; 1(5):769-71.
3. *Numerical Recipes* for subroutines MNBRAK and FMAXIM.

## Procedure Outline

1. **Principal Component Analysis (PCA)**  
   Conduct PCA on the Intrinsic Reaction Coordinate (IRC) in Cartesian space, representing atomic motions linking various stationary points on the PES. Use [PathReducer](https://github.com/share1992/PathReducer) for efficient PCA on IRC data.

2. **Wavelet Transform**  
   Apply the wavelet transform to the reduced data from PCA to analyze time-dependent frequencies and visualize atomic motions over time.

3. **Spline Fitting (Optional)**  
   For further analysis, you can perform spline fitting on the wavelet transform results using `spline.py`. This requires Python with `scipy` and `numpy` installed.

---

## Setup and Requirements

### Required Files
1. **wavelet.f** - Program to calculate the wavelet transform.
2. **input.dat** - Contains time vs. coordinates along the trajectory.
3. **norm.gnu** - Gnuplot settings file to plot the 3D wavelet figure.

### Output Files
Running the program will generate these output files:
1. **wavelet.exe** - Executable file created from `wavelet.f`.
2. **fort.66** - Contains the non-normalized 3D data.
3. **fort.8** - Contains the normalized 3D data in a gnuplot-plottable format.
4. **max.out** - Contains the maximum amplitude frequency over time.
5. **3D.pdf** - 3D wavelet plot.

### Adjustable Parameters in `wavelet.f`
1. `N` - Number of trajectory points for the wavelet calculation, equivalent to the total number of lines in `input.dat` (update `N` in four places).
2. `BS` - Every BS-th line from `input.dat` is used for the wavelet (default: `BS=10`).
3. `BS1` - Number of time points computed, given by `N/BS` (update in the declaration section).
4. `FRR` - Number of frequency points per time point. The total frequency range is from `A0` to `A0 + 10*FRR` (default: `FRR=500`).

   > **Note**: Defaults for `BS` and `FRR` are typically adequate. The `BS` value is associated with frequency: low-frequency spectral information should have a larger wavelet time interval, and vice versa.

### Adjustable Parameters in `input.dat`
1. **Flag**: The first line is a flag for printing 3D data. If `dim3=1`, 3D data is printed; otherwise, it is not.
2. **Minimum Frequency (`A0`)**: The second line is `A0`, the minimum frequency for starting calculations. Ensure that low-frequency artifacts have diminished before `A0`. To check, conduct wavelet transforms at selective time points.

---

## Installation and Usage

### Software Requirements
- Fortran compiler (e.g., `gfortran`)
- Gnuplot (for generating 3D plots)
- Python with `scipy` and `numpy` (for spline fitting)

### Steps to Run
1. Compile the Fortran code:
    ```bash
    gfortran wavelet.f -o wavelet.exe
    ```
2. Run the program with:
    ```bash
    ./wavelet.exe < input.dat > max.out
    ```
3. Plot the data with gnuplot:
    ```bash
    gnuplot norm.gnu
    ```
4. For visualizing `max.out` in gnuplot:
    ```gnuplot
    plot 'max.out' using 1:2 with lines
    ```

---

## Contributors
- **Manikandan Paranjothy** – Primary contributor
- **Richa Rashmi** – Contributed modifications
- **Aniruddha Seal** – Contributed modifications

---

## Citations
If you use this code in your research, please cite the following works:

```bibtex
@article{rashmi2023z,
  title={E- Z Isomerization in Guanidine: Second-order Saddle Dynamics, Non-statisticality, and Time-frequency Analysis},
  author={Rashmi, Richa and Yadav, Pankaj Kumar and Seal, Aniruddha and Paranjothy, Manikandan and Lourderaj, Upakarasamy},
  journal={ChemPhysChem},
  volume={24},
  number={2},
  pages={e202200640},
  year={2023},
  publisher={Wiley Online Library}
}
```

```bibtex
@article{manikandan2014dynamical,
  title={Dynamical traps lead to the slowing down of intramolecular vibrational energy flow},
  author={Manikandan, Paranjothy and Keshavamurthy, Srihari},
  journal={Proceedings of the National Academy of Sciences},
  volume={111},
  number={40},
  pages={14354--14359},
  year={2014},
  publisher={National Acad Sciences}
}
```

```bibtex
@article{hare2019low,
  title={Low dimensional representations along intrinsic reaction coordinates and molecular dynamics trajectories using interatomic distance matrices},
  author={Hare, Stephanie R and Bratholm, Lars A and Glowacki, David R and Carpenter, Barry K},
  journal={Chemical science},
  volume={10},
  number={43},
  pages={9954--9968},
  year={2019},
  publisher={Royal Society of Chemistry}
}
```
