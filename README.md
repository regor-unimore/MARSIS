# A Predict-then-Optimize Approach for the Research of Underground Water on Mars
This repository contains the companion code and instances for the paper "A Predict-then-Optimize Approach for the Research of Underground Water on Mars" by Benedetta Ferrari <sup>[1]</sup>, Maxence Delorme <sup>[2]</sup>, Manuel Iori <sup>[1]</sup>, Marco Lippi <sup>[3]</sup>, and Roberto Orosei <sup>[4]</sup>.

---

**Instances** contains a folder for each of our test instances, named as M_G_L_NC, where G is the granularity of the discretization, L the length of the observation, and NC the number of quality classes, as defined in the article. In each folder there are four .txt files, each corresponding to a different component of the instance:
<ol>
    <li>list of the PIs and related features. Legend of columns:
        <ul>
            <li>Id</li>
            <li>x coordinate (m)</li>
            <li>y coordinate (m)</li>
            <li>PIAV: 1 if PI was already covered in the past; 0 otherwise</li>
            <li>initial quality (dB)</li>
            <li>initial quality class</li>
            <li>s_p: 1 if PI belongs to an AOI; 0 otherwise</li>
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
    <li> <i> ephemeris_time.txt </i> : sampled times for future observations; each observation starts at a sampled time and last for at most L seconds</li>
    <li> <i> predictions.txt </i> : machine learning predictions for each sampled time</li>
</ul>

---

[1] Department of Science and Methods for Engineering, University of Modena and Reggio Emilia, Via Giovanni Amendola 2, Reggio Emilia 42122, Italy. <br>
[2] Department of Econometrics and Operations Research, Tilburg University, 5037 AB Tilburg, The Netherlands. <br>
[3] Department of Information Engineering, University of Florence, Via di Santa Marta, 3, Florence, 50139, Italy. <br>
[4] Institute of Radioastronomy, Italian National Institute of Astrophysics, Via Piero Gobetti 101, Bologna 40129, Italy. <br>
