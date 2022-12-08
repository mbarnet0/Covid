# Covid
Code for PNAS Paper "Epidemic Responses Under Uncertainty" by Michael Barnett, Greg Buchak, and Constantine Yannelis.

To derive the numerical results in the paper, there are three sets of code that need to be run:
1. Outside-the-Model Uncertainty Code
2. Inside-the-Model Uncertainty Code
3. Detection Error Probability Code

This repository contains a folder for each setting. I outline the code provided, as well as details about to run the code and in what order the code needs to be run, in order to allow users to replicate the results provided in the paper.

## Outside-the-Model Uncertainty Code
There are three separate codes to run for this setting:
1. BBY_PNAS_COVID19_OutsideUncertaintyCases.m
2. BBY_PNAS_COVID19_OutsideUncertaintySmartCases.m
3. BBY_PNAS_COVID19_OutsideUncertainty_Sims.m

## Inside-the-Model Uncertainty Code
There are three separate codes to run for this setting:
1. BBY_PNAS_COVID19_InsideUncertaintyCases.m
2. BBY_PNAS_COVID19_InsideUncertaintySmartCases.m
3. BBY_PNAS_COVID19_InsideUncertainty_Sims.m

## Detection Error Probability Code
There are two separate codes to run for this setting:
1. BBY_PNAS_COVID19_DEPs.m
2. BBY_PNAS_COVID19_DEPs_Comp.m

## Replicating the Results

In order to replicate the results, you will need to run each of these codes for each of the cases specified in the paper. There are 16 outside-the-model uncertainty cases, 8 inside-the-model uncertainty cases, and 4 sets of DEPs for each of the ``uncertainty averse'' inside-the-model uncertainty cases. The code in each setting needs to be run in the order listed above, where the ``SmartCases'' codes are optional and can be used to derive solutions with alternative parameter values using existing solutions as initial guesses for the value function.

Each model solution case takes various days to run. We ran the jobs for each setting simultaneously on the ASU Sol and Agave high-performance-computation servers. Moreover, the DEP computation is such that it can be significantly parallelized, which we exploited by splitting the daily observations for the model solution simulations into simultaneous batch jobs and parallelized each batch job.
