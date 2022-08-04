# MVF_Version2
MVF (modal based vector fitting) is an analytical method for coupling matrix extraction of microwave filters. Version 2 modifies the code to make it faster and more reasonable. 

Please run the code in MATLAB newer than R2018b.
1. Run MVF_main.m;
2. Import the simulated/measured S-parameters as data file in .s2p format;
3. Output the extracted coupling matrix representing the equivalent circuit.

References：
1. Circuit Model Extraction of Parallel-Connected Dual-Passband Coupled-Resonator Filters (for fit Y-parameters with rational functions as pole-residue form)
2. Phase De-Embedding of Narrowband Coupled-Resonator Networks by Vector Fitting (for phase de-embedding of the S-parameters)
