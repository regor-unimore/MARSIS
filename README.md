This repository contains the tested instances and code for all algorithms discussed in the paper "A Predict-then-Optimize Approach for the 
Research of Underground Water on Mars" by Benedetta Ferrari, Maxence Delorme, Manuel Iori, Marco Lippi and Roberto Orosei. 

The material is divided into "Machine Learning" and "Optimization" data/code, to better distinguish the two part of the proposed predict-then-optimize approach. 
All algorithms are coded in Python and the ILP model require the commercial solver Gurobi (we used version 11.0.1). In the followinf the content of each folder is described:

*Machine Learning (ML)*

**Data** contains the link to access to the historical dataset and the list of orbit to exclude during training as they were considered compromised by solar events. It also contains the dataset for "future" planning, used to generate the predictions and the optimization istances (MARSIS_dataset_2023_2025).

**Code** contains the algorithms setting used to train the five regression models, based on scikit-learn package and keras (for the neural network)

*Optimization (Opt)*

**Instances** contains a folder for each of our test instances, named as M_G_L_NC, where G is the granularity of the discretization, L the length of the observation, and NC the number of quality classes, as defined in the article. In each folder there are four .txt files, each corresponding to a different component of the instance:
<ol>
    <li>list of the PIs and related features. Legend of columns:
        <ul>
            <li>Id of the PI</li>
            <li>x<sub>p</sub>: x coordinate (m)</li>
            <li>y<sub>p</sub>: y coordinate (m)</li>
            <li>cov<sub>p</sub>: 1 if PI was already covered in the past; 0 otherwise</li>
            <li>q<sub>p</sub>: initial quality (dB)</li>
            <li>c<sub>p</sub>: initial quality class</li>
            <li>w<sub>p</sub>: 1 if PI belongs to an AOI; 0 otherwise</li>
        </ul>
    </li>
    <li>list of PIs covered by each observation (observation = row)</li>
    <li>list of continuous quality associated to the PIs in each observation (machine learning predictions)</li>
    <li>list of quality classes associated to the PIs in each observation</li>
</ol>

**General** contains files that can be used in combination with each instance. Specifically:
<ul>
    <li> <i> SouthPole_contour.dat </i> : list of (x,y) defining the contour of the South Pole</li>
    <li> <i> AOI*_countour.dat </i> : lists of (x,y) defining the four areas of interest</li>
    <li> <i> day.txt </i> : date of each observation (observation = row)</li>
    <li> <i> orbit.txt </i> : orbit id of each observation (observation = row)</li>
    <li> <i> ephemeris_time.dat </i> : sampled times for future observations; each observation starts at a sampled time and last for at most L seconds</li>
    <li> <i> predictions.txt </i> : machine learning predictions for each sampled time</li>
    <li> <i> x.dat </i> : x coordinate for each sampled time</li>
    <li> <i> y.dat </i> : y coordinate for each sampled time</li>
</ul>

**Code** contains the python code to solve the ILP model with Gurobi solver and a visualization template for the model solution. The codes take in input the name of the instance to be solved (or visualized), and read all the necessary input files. Note that in the ILP N, W, and &beta; are set to their deafult values, and can be changed directly in the code.

**Computational Results**: computational results for all tested algorithms




