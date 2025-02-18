# README

## MATLAB Scripts for Generating Figures

This repository contains MATLAB scripts used to generate figures in the manuscript. Below is a description of each script and the corresponding figures they produce.

### Scripts and Their Corresponding Figures

#### **1. Runner1.m**
Used to generate Figures 2 and 3.

- **Figure 2**:
  - Transmission functions for both models.
  - Fitness function \( \phi_1 \) for Model 1.
  - Transmission functions: \( \alpha(x) = 1.5\sech^2(10x - 0.2) \) and \( \beta(x) = \sech^2(10x - 2) \), based on death rate \( x \).
  - The fitness function \( \phi_1 \) attains its maximum at \( \mu_1^* = 0.17435 \), representing the Evolutionarily Stable Strategy (ESS) value for Model 1.

- **Figure 3**:
  - Fitness function \( \phi_2 \) for the model incorporating asymptomatic transmission.
  - (a) 2D domain of \( \phi_2(\mu_A, \mu_S) \).
  - (b) \( \phi_2 \) plotted with \( \mu_A \) as the evolutionary variable, with different \( \mu_S \) values. The solid red line represents the critical \( \mu_S \) value that maximizes fitness.
  - (c) \( \phi_2 \) plotted with \( \mu_S \) as the evolutionary variable, with different \( \mu_A \) values. The solid red line represents the critical \( \mu_A \) value that maximizes fitness.

#### **2. sensityvityNu.m**
Used to generate Figure 4.

- **Figure 4**:
  - Effect of mean infectious period (\( \nu^{-1} \)) on ESS and maximum fitness values for both models.
  - (a) ESS values \( \mu_1^* \), \( \mu_A^* \), and \( \mu_S^* \) as a function of \( \nu^{-1} \).
  - (b) Maximum fitness values for both models, showing greater stability for Model 2.
  - (c) and (d) Detailed views of \( \mu_1^* \) and \( \mu_S^* \) with respect to \( \nu^{-1} \).
  - (e) Relationship between \( \mu_A^* \) and \( \nu^{-1} \) for Model 2.

#### **3. sensitivityP.m**
Used to generate Figures 5a-5c.

#### **4. SEnsityVityOmega.m**
Used to generate Figures 5d-5f.

- **Figure 5**:
  - Effect of the expected time in the asymptomatic state (\( \omega^{-1} \)) and fraction of individuals following the "mild" recovery pathway (\( p \)) on ESS values and maximum fitness.
  - (a) ESS values \( \mu_A^* \) and \( \mu_S^* \) as functions of \( \omega^{-1} \).
  - (b) Detailed behavior of \( \mu_A^* \) with respect to \( \omega^{-1} \).
  - (c) Maximum fitness values for Model 2 as \( \omega^{-1} \) increases, surpassing Model 1 at \( \omega^{-1} = 3 \).
  - (d) ESS values \( \mu_A^* \) and \( \mu_S^* \) as functions of \( p \).
  - (e) Detailed behavior of \( \mu_A^* \) with respect to \( p \).
  - (f) Maximum fitness values for Model 2 decrease with increasing \( p \), converging to Model 1 as \( p \to 1 \).

### Usage Instructions
1. Ensure you have MATLAB installed.
2. Run the scripts in MATLAB to generate the respective figures.
3. Adjust parameter values within scripts if needed to explore different scenarios.

For any issues or inquiries, please refer to the manuscript or reach out to the repository maintainer.
