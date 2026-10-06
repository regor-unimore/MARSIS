This archive is distributed in association with the journal Operations Research under the MIT license. 

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper "A Predict-then-Optimize Approach for the 
Research of Underground Water on Mars" by B. Ferrari, M. Delorme, M. Iori, M. Lippi and R. Orosei.

Important: This code is being developed on an on-going basis at https://github.com/](https://github.com/regor-unimore/MARSIS. Please go there if you would like to get a more recent version or would like support.

***Cite***
--

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.
[to be completed]

***Description***
--

The material within the data and script folders is divided into "Machine Learning" and "Optimization", to better distinguish the two parts of the proposed predict-then-optimize approach. In the following, the content of each folder is described:

-------

*Machine Learning (ML)*

**Data** contains the link to access to the historical dataset and the list of orbits to exclude during training, as they were compromised by solar events. It also contains the dataset for "future" planning, used to generate the predictions and the optimization istances (MARSIS_dataset_2023_2025).

**Code** contains the algorithms setting used to train the five regression models based on scikit-learn package (**Test_KFold.py**) and the keras-based neural network (**Test_KFold_DNN_keras.py**).

-------

*Optimization (Opt)*

**Instances** contains a folder for each of our test instances, named as M_G_L_NC, where G is the granularity of the discretization, L the length of the observation, and NC the number of quality classes, as defined in the article. Each folder contains four .txt files, each corresponding to a different component of the instance:
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

**Code** contains the python code to solve all the optimization algorithms tested within the research. The codes take in input the name of the instance to be solved and read all the necessary input files. In each code, the solution parameters are set to their default values and can be changed directly in the code. Moreover, within the matheuristics codes (type R, C, and G) it is possible to chose which variant to execute. A further code is included to help visualize the solution.

***Requirements and Execution***
--

All algorithms are coded in Python and the ILP model requires the commercial solver Gurobi (we used version 11.0.1). The experiments were done on a single thread running on an Intel(R) Xeon(R) CPU E5-2640 v3 at 2.60 GHz under Windows 10 Pro with 64 GB RAM. 

--> To train the regression models, first download the historical dataset and the list of orbits to be excluded. Place both files in the same directory as the source code. Once the data are correctly set up, the training scripts can be executed!

--> To run the optimization algorithms, download the instance folders and the general data files, and place them in the same directory as the source code. Once the data are properly set up, the code can be executed. Note that, if necessary, default parameters (e.g., N and W) have to be manually modified within the code before execution.

***Execution order and computational cost***
--

The heuristic H2 must be run before the Lagrangian relaxation and the matheuristics, since they use the H2 solution as initial solution (by default, they read it from `output/H2_<instance>.txt`):

```
setup_data.py  -->  H1.py, H2.py, ILP Model.py, Benders.py         (independent of each other)
                    H2.py  -->  Lagrangian.py, Matheuristic_R.py, Matheuristic_C.py, Matheuristic_G.py
```

The table below reports the computational effort required to reproduce each sheet of `results/Computational Results.xlsx`, computed from the times reported in the workbook (single thread, hardware described above). All runs are single-threaded and independent of each other, so they can be executed in parallel.

| Experiment (workbook sheet) | Runs | Time limits | Total time |
|---|---|---|---|
| ILP model (*ILP*) | 40 | 3600 s per run | about 24 h (15 runs reach the time limit; 7 runs on G = 01 instances report no time) |
| Heuristic H1 (*Heuristics*) | 40 | none | about 15 min (at most 2 min per run) |
| Heuristic H2 (*Heuristics*) | 40 | none | about 7.6 h (at most 2 h per run, on M_01_268_4) |
| Lagrangian relaxation (*Lagrangian*) | 40 | none | not reported in the workbook (1-2 min per run on G = 10 instances) |
| Matheuristics R, C and G (*Matheuristic all*) | 5600 (10 or 20 repetitions per configuration) | 3600 s per repetition, 360 s per Gurobi call | about 3900 h |

A single execution of a matheuristic script with the default settings performs 20 repetitions (`--n-iter`) and can therefore take up to about 20 hours.
