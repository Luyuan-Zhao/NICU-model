# Optimising Intervention Strategies for ESBL-Enterobacterales in a NICU: An Agent-Based Model


This repository contains the source code for the agent-based model (ABM) developed for the dissertation titled: **"Optimising Intervention Strategies for Extended-Spectrum Beta-Lactamase-Producing *Enterobacterales* in a Tanzanian Neonatal Intensive Care Unit: An Agent-Based Modelling Approach"**.

Author: **Luyuan Zhao** 
Supervisor: **Professor Louise Matthews** 

---

## Overview

In Neonatal Intensive Care Units (NICUs), particularly in low- and middle-income countries (LMICs), the spread of Extended-Spectrum Beta-Lactamase (ESBL)-producing *Enterobacterales* is a major public health threat, contributing significantly to neonatal morbidity and mortality.

This project introduces an ABM designed to simulate the complex transmission dynamics of these pathogens within a resource-limited, 10-bed NICU in Mwanza, Tanzania. The model quantitatively evaluates the effectiveness of various Infection Prevention and Control (IPC) measures and antimicrobial stewardship strategies, both individually and in combination.

The primary goal is to provide a data-driven tool to help identify the most effective and resource-efficient intervention bundles for controlling antimicrobial resistance (AMR) in high-burden settings.

## Key Findings from the Study

Our simulations revealed several key insights into AMR dynamics and control:

1.  **Cross-Transmission is the Main Amplifier**: While antibiotic selection pressure is the catalyst for resistance, contact-mediated cross-transmission (e.g., via healthcare workers or contaminated equipment) is the primary engine driving the spread of resistant strains. Transmission events vastly outnumbered within-host selection events.
2.  **Hierarchy of Interventions**: Not all interventions are equally effective.
    * **Most Impactful Single Actions**: Reducing empiric antibiotic use (Strategy A) and enhancing the cleaning of shared medical equipment (Strategy D) were the most effective standalone strategies.
    * **The Hand Hygiene Paradox**: Improving hand hygiene alone had a minimal impact in the model, suggesting its effectiveness is limited when other redundant transmission pathways (like contaminated surfaces) are not addressed.
3.  **Diminishing Returns in Bundled Strategies**: A comprehensive bundle of all interventions (Strategy H) was the most effective overall, reducing cumulative resistant infections by ~54%. However, it exhibited a slight antagonistic effect (diminishing returns), meaning its total benefit was less than the sum of its individual parts.
4.  **A "Lean, High-Impact Bundle" is Proposed**: The findings suggest that a pragmatic and resource-efficient "lean bundle"—prioritising **robust antibiotic stewardship** and **targeted cleaning of high-touch medical equipment**—may yield the greatest health benefits at a lower implementation cost.

## Model Features

The ABM simulates the interactions between three types of agents: **Neonates**, **Healthcare Workers (HCWs)**, and **Environmental Reservoirs**.

* **Detailed Neonate States**: Neonates transition between six states: Susceptible, Colonised (sensitive/resistant), Infected (sensitive/resistant), and Recovered.
* **Multiple Transmission Pathways**: The model integrates direct contact (HCW-to-neonate) and indirect contact (via contaminated medical devices and general environmental surfaces).
* **Antibiotic Pressure**: Simulates how empiric antibiotic use creates selection pressure, potentially causing a within-host shift from sensitive to resistant strains.
* **Realistic NICU Dynamics**: Incorporates stochastic patient admission/discharge, a fixed HCW-to-neonate ratio, and is parameterized using real-world longitudinal surveillance data from a Tanzanian NICU.
* **Factorial Simulation Framework**: Allows for testing and comparing 8 different single and bundled intervention strategies across various baseline transmission intensities.


---

## System Requirements & Installation

The simulations and analysis require both **Julia** and **R**.

### 1. Julia (for running the ABM simulations)

* **Language**: Julia v1.10.5 .
* **Packages**:
    * `Random`
    * `StatsBase` 
    * `DataFrames` 
    * `Plots` 

You can install these packages in the Julia REPL:
```julia
using Pkg
Pkg.add("Random")
Pkg.add("StatsBase")
Pkg.add("DataFrames")
Pkg.add("Plots")
```

### 2. R (for statistical analysis and plotting)

* **Language**: R v4.4.1 
* **Packages**:
    * `ggplot2` 
    * `dplyr` 
    * `tidyr` 
    * `patchwork` 
    * `epitools` 

You can install these packages in the R console:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "epitools"))
```

## How to Use

All simulation and analysis scripts are available in this repository to ensure full reproducibility.

1.  **Run the ABM Simulations**:
    * Navigate to the `/scripts` (or relevant) folder.
    * Execute the main Julia script to run the simulations. The model supports running single scenarios or the full 24-scenario factorial analysis (8 strategies x 3 transmission intensities).
    * *(1. Main0713-V1.14.jl,2. ScenarioSimulations_0713.jl,3. SensitivityAnalysis_0713 .jl are the main function to run the ABM in Julia, then the .csv files can be generated automatically as listed here (e.g. all_scenarios_full.csv, sensitivity_X.csv))*

2.  **Analyse the Results**:
    * The simulation scripts will generate output files (e.g., CSVs) in the `/results` (or relevant) folder.
    * Execute the R scripts to process these output files, generate summary statistics, and create the plots (like Figures 3, 4, and 5) from the paper.
    * *(All stats used in R are consistent with its names, which means you can directly run the codes)*


