# Optimising Intervention Strategies for ESBL-Enterobacterales in a NICU: An Agent-Based Model


[cite_start]This repository contains the source code for the agent-based model (ABM) developed for the dissertation titled: **"Optimising Intervention Strategies for Extended-Spectrum Beta-Lactamase-Producing *Enterobacterales* in a Tanzanian Neonatal Intensive Care Unit: An Agent-Based Modelling Approach"**[cite: 3, 8].

[cite_start]Author: **Luyuan Zhao** [cite: 4]
[cite_start]Supervisor: **Professor Louise Matthews** [cite: 6]

---

## Overview

[cite_start]In Neonatal Intensive Care Units (NICUs), particularly in low- and middle-income countries (LMICs), the spread of Extended-Spectrum Beta-Lactamase (ESBL)-producing *Enterobacterales* is a major public health threat, contributing significantly to neonatal morbidity and mortality[cite: 10].

[cite_start]This project introduces an ABM designed to simulate the complex transmission dynamics of these pathogens within a resource-limited, 10-bed NICU in Mwanza, Tanzania[cite: 11, 12]. [cite_start]The model quantitatively evaluates the effectiveness of various Infection Prevention and Control (IPC) measures and antimicrobial stewardship strategies, both individually and in combination[cite: 11].

[cite_start]The primary goal is to provide a data-driven tool to help identify the most effective and resource-efficient intervention bundles for controlling antimicrobial resistance (AMR) in high-burden settings[cite: 20, 21].

## Key Findings from the Study

Our simulations revealed several key insights into AMR dynamics and control:

1.  [cite_start]**Cross-Transmission is the Main Amplifier**: While antibiotic selection pressure is the catalyst for resistance, contact-mediated cross-transmission (e.g., via healthcare workers or contaminated equipment) is the primary engine driving the spread of resistant strains[cite: 15, 259, 349]. [cite_start]Transmission events vastly outnumbered within-host selection events[cite: 258, 259].
2.  **Hierarchy of Interventions**: Not all interventions are equally effective.
    * [cite_start]**Most Impactful Single Actions**: Reducing empiric antibiotic use (Strategy A) and enhancing the cleaning of shared medical equipment (Strategy D) were the most effective standalone strategies[cite: 296, 298].
    * [cite_start]**The Hand Hygiene Paradox**: Improving hand hygiene alone had a minimal impact in the model, suggesting its effectiveness is limited when other redundant transmission pathways (like contaminated surfaces) are not addressed[cite: 368, 452].
3.  [cite_start]**Diminishing Returns in Bundled Strategies**: A comprehensive bundle of all interventions (Strategy H) was the most effective overall, reducing cumulative resistant infections by ~54%[cite: 16, 293]. [cite_start]However, it exhibited a slight antagonistic effect (diminishing returns), meaning its total benefit was less than the sum of its individual parts[cite: 18, 364].
4.  [cite_start]**A "Lean, High-Impact Bundle" is Proposed**: The findings suggest that a pragmatic and resource-efficient "lean bundle"—prioritising **robust antibiotic stewardship** and **targeted cleaning of high-touch medical equipment**—may yield the greatest health benefits at a lower implementation cost[cite: 20, 461, 462].

## Model Features

[cite_start]The ABM simulates the interactions between three types of agents: **Neonates**, **Healthcare Workers (HCWs)**, and **Environmental Reservoirs**[cite: 13, 117].

* [cite_start]**Detailed Neonate States**: Neonates transition between six states: Susceptible, Colonised (sensitive/resistant), Infected (sensitive/resistant), and Recovered[cite: 125].
* [cite_start]**Multiple Transmission Pathways**: The model integrates direct contact (HCW-to-neonate) and indirect contact (via contaminated medical devices and general environmental surfaces)[cite: 13, 115, 132].
* [cite_start]**Antibiotic Pressure**: Simulates how empiric antibiotic use creates selection pressure, potentially causing a within-host shift from sensitive to resistant strains[cite: 13, 130].
* [cite_start]**Realistic NICU Dynamics**: Incorporates stochastic patient admission/discharge, a fixed HCW-to-neonate ratio, and is parameterized using real-world longitudinal surveillance data from a Tanzanian NICU[cite: 12, 103, 141, 142].
* [cite_start]**Factorial Simulation Framework**: Allows for testing and comparing 8 different single and bundled intervention strategies across various baseline transmission intensities[cite: 14, 158, 159].


---

## System Requirements & Installation

The simulations and analysis require both **Julia** and **R**.

### 1. Julia (for running the ABM simulations)

* [cite_start]**Language**: Julia v1.10.5 [cite: 150]
* **Packages**:
    * `Random`
    * [cite_start]`StatsBase` [cite: 151]
    * [cite_start]`DataFrames` [cite: 151]
    * [cite_start]`Plots` [cite: 151]

You can install these packages in the Julia REPL:
```julia
using Pkg
Pkg.add("Random")
Pkg.add("StatsBase")
Pkg.add("DataFrames")
Pkg.add("Plots")
```

### 2. R (for statistical analysis and plotting)

* [cite_start]**Language**: R v4.4.1 [cite: 186]
* **Packages**:
    * [cite_start]`ggplot2` [cite: 194]
    * [cite_start]`dplyr` [cite: 194]
    * [cite_start]`tidyr` [cite: 194]
    * [cite_start]`patchwork` [cite: 194]
    * [cite_start]`epitools` [cite: 193]

You can install these packages in the R console:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "epitools"))
```

## How to Use

[cite_start]All simulation and analysis scripts are available in this repository to ensure full reproducibility[cite: 151, 195].

1.  **Run the ABM Simulations**:
    * Navigate to the `/scripts` (or relevant) folder.
    * Execute the main Julia script to run the simulations. [cite_start]The model supports running single scenarios or the full 24-scenario factorial analysis (8 strategies x 3 transmission intensities)[cite: 159].
    * *(1. Main0713-V1.14.jl,2. ScenarioSimulations_0713.jl,3. SensitivityAnalysis_0713 .jl are the main function to run the ABM in Julia, then the .csv files can be generated automatically as listed here (e.g. all_scenarios_full.csv, sensitivity_X.csv))*

2.  **Analyse the Results**:
    * [cite_start]The simulation scripts will generate output files (e.g., CSVs) in the `/results` (or relevant) folder[cite: 152].
    * Execute the R scripts to process these output files, generate summary statistics, and create the plots (like Figures 3, 4, and 5) from the paper.
    * *(You may want to add a specific command here, e.g., `Rscript analyze_results.R`)*


