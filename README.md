# Usage: 
<pre>
gdilp.py [-h] [--intervallimit [INTERVAL_LIMIT]]
         [--max_converted_spec [MAX_CONVERTED_SPEC]]
         [--time_limit [TIME_LIMIT]]
         input_file
 </pre> 

### positional arguments:

`input_file`  &ensp;&ensp; sets the name of the file with the input data;<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Default input is a file in RME format. The first line consist of the number of gene trees n. <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; The next n lines represent gene trees. Then, the following line consist of the species tree.<br>

### optional arguments:

`-h, --help`  &ensp;&ensp;shows help message and exit<br>

`--intervallimit [INTERVAL_LIMIT]` <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; sets the maximal length of intervals (parameter _k_); <br> 
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Parameter _k_ is set to INTERVAL_LIMIT value.<br> 
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;  Default _k = 0_, and no intervals are shortened. It is equivalent to _&lambda; = &infin;_ used in the publication.<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; If _k > 0_, then all intervals with length _y > k_ are shortened to the size of _k_. <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Note _&lambda;_ from the definition in the publication (edge interval length) <br> 
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; equals _k-1_ (as parameter _k_ denote node length of an interval).<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; In example for a gene tree node _g_ the interval equals _&lang;M(g),b&rang;_, where _M(g)_ is a lca-mapping of _g_<br> 
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;  in the species tree, and _b_ is a node in species tree.<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; If _|&lang;M(g),b&rang;| = y > k_ then we assign new interval _&lang;M(g),c&rang;_ such that _|&lang;M(g),c&rang;| = k_ and _c &isin; &lang;M(g),b&rang;_.<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;  In other words, nodes are removed from 'b' side of the interval. <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; In summary, for every gene tree node _g_ we limit the distance (in the species tree) between <br> 
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; its mapping in a potential solution and its lca-mapping  not to exceed _k_.<br>

`--max_converted_spec [MAX_CONVERTED_SPEC]` <br> 
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; sets the maximal number of speciations that can be converted into duplications (parameter _&sigma;_).<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Parameter _&sigma;_ is set to MAX_CONVERTED_SPEC value.<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Default value of the parameter is _&sigma; = &infin;_. To manually set _&sigma; = &infin;_, write -1. <br>

`--time_limit [TIME_LIMIT]` sets the maximal computation time in seconds.<br>



# Installation:

### Step 1 - Download sources

Download zipped folder of this repository. Unzip to project source folder (like `genomicduplicationilp-master/`).<br>

### Step 2 A - Run on a local machine

**1.** **Install gurobi solver.** <br>
**-** see https://www.gurobi.com/downloads/gurobi-software/ <br>
**-** note, it is possible to obtain an academic license <br>
**-** gurobipy package can be installed by `python3.9 -m pip install gurobipy`<br>
**-** gurobi licence can be installed by `grbgetkey <licence number>` <br>
**-** gurobipy version 9.1.1 was used

**2.** **Install other standard python packages if they are missing.** <br>
**-** See `requiremets.txt` for other requirements info. <br>
**- Requirements** <br>
Default modules: &emsp; &emsp;&emsp;&ensp;&ensp;timeit,  argparse (Python >= 3.2), io<br>
Required modules: &emsp; &emsp; &ensp; Bio, gurobipy, pandas<br>
Optional:&emsp; &emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&ensp;pytest<br>
Sample commands to install required modules: <br>
**- -** `python3.9 -m pip install Biopython` <br>
**- -** `python3.9 -m pip install pandas` <br>

**3.** **Run** <br>
**-** Sample command`python gdilp.py input_file`

### Step 2B - Run on a remote server

The installation of required packages (like gurobipy) may be possible only with root rights. <br> 
Therefore, I suggest to create a python environment as described below.<br>
The following steps could be also used on local machines.

### Step 2B.1 - First run (with the creation of an environment)
**2B.1.1. Create python environment**<br>

(1) `python3 -m venv gurobi-env` &emsp;&emsp; &emsp;&emsp;&emsp;&ensp;&emsp; Create virtual environment at project source folder.<br>
(2) `source ./gurobi-env/bin/activate` &emsp;&emsp;&emsp;&ensp; Activate the environment.<br>

**2B.1.2. gurobipy installation:**<br>
(3) `cd /opt/gurobi901/linux64/` &emsp;&emsp;&ensp;&emsp;&emsp;&emsp;&emsp;&ensp; Go to the gurobi folder.<br>
(4) `python3 setup.py build -b PATH install` &ensp; Gurobi installation to a selected PATH. Useful on remote servers. <br>
_Note: Running the above command not in virtual environment (without steps 1-2) is not enough, because it instals <br> 
gurobipy into the python folder causing another permission denied error (not in /opt/gurobi... but in python folder instead)_<br>
(5) `cd genomicduplicationilp-master/` &emsp;&emsp;&emsp;&ensp; Go to the project source folder.<br>

**2B.1.3. Install other standard python packages if they are missing:**<br>
(6) `pip3 install biopython ` &emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&emsp; &emsp;&emsp;&ensp; Biopython installation.<br>
(7) `pip3 install pandas` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;&ensp; Pandas installation.<br>

**Now the environment is ready to use. Use commands (ii)-(iv). Then multiple runs like (v). See below.**<br>
**To close virtual environment use (vi).**

### Step 2B.2 - Run from environment

**2B.2.1 Activate the environment.** <br>
( i ) `source ./gurobi-env/bin/activate` &emsp;&emsp; &emsp;&emsp;&emsp;&emsp;<br>
**2B.2.2 Configure the environment.** <br>
(ii )  `export GUROBI_HOME="/opt/gurobi901/linux64"` Now 3 magic configuration lines.<br>
(iii)  `export PATH="${PATH}:${GUROBI_HOME}/bin"`<br>
(iv) `export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"`<br>
**2B.2.3 Now multiple runs are possible, in example:** <br>
(v ) `python3 gdilp.py tests/guigo.txt` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;<br>
**2B.2.4 To close the environment:**<br>
(vi)  `deactivate`&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp;&emsp;<br>

Note:<br>
_Normally, we only need to activate the environment. And than run as many computations (v) as we desire._<br>
_However, if we log out from the server and then log in again, then we have to run 3 magic lines again.<br>_

## Project files

**Simulated data**<br>
The original SimPhy results are in `wgd-simulations` folder.<br>
`sim_utils.py` was used to extract data to input files stored as `tests\sim[i].txt`<br>
Note `[i]`file corresponds to _S<sub>i</sub>_ tree (see Figure 1)<br>

**Input data**<br>
Input data is stored in `tests\` folder. In our study we use TreeFam dataset - file `treefam.txt`.<br>
The folder contains also sample tests and `test_gd_ilp.py` to perform unittest.<br>

**Output data**<br>
The experimental section describes the evaluation of the results stored at `results/out_name_k_`&sigma;`.csv`<br>
Where:<br>
- name - sim1, ..., sim5, treefam - is the name of the corresponding input dataset<br>
- k - the value of INTERVAL_LIMIT parameter used (default or given by user)<br>
- &sigma; -the value of MAX_CONVERTED_SPEC parameter used (default or given by user)<br>

## Input / Output formats:

### Input file <br>
Default input is a file in the following format. The first line consist of the number of gene trees n. <br>
The next n lines represent gene trees. Then, the following line consist of the species tree.<br>

### Output <br>

**1. Console output** <br>
-- First we display chosen parameters <br>
-- Next, ILP model details are shown<br>
-- Last lines consist of:<br>
---- dictionary of the gene tree support for duplication events (see below)<br>
---- Total MEscore value <br>
---- Time of the computation <br>
_Note, to save above info into save_info_file, run command with `>save_info_file`<br>
in example`python3 gdilp.py tests/guigo.txt >save_info_file`_ <br>

**2. Detailed output** is stored in`results/out_name_k_`&sigma;`.csv`, where <br>
k and &sigma; are parameter values (see above), and _name_ is an input file name. <br>
-- Each row in csv file corresponds to a node _s_ in a species tree. <br> 
-- Columns in the csv file represent:<br>
---- Node Cluster - names of all species, which are descendants of _s_; separated by semicolons<br> 
&emsp; **!!! Note, do not choose semicolon as separator while openning csv file, choose only comma** <br> 
---- ME score - that corresponds to node _s_; <br>
&emsp; _ME score denotes number of comparable duplication clusters mapped to s._ <br>
---- Lets assume ME score for _s_ equals 2. Then we have 2 non-empty pairs of columns:<br>
&emsp;&emsp; Level 0 support - number of gene trees supporting one duplication event at _s_<br>
&emsp;&emsp; Level 0 trees - set of id numbers of gene trees supporting one duplication event at _s_<br>
&emsp;&emsp; Level 1 support - number of gene trees supporting two comparable duplication events at _s_<br>
&emsp;&emsp; Level 1 trees - set of id numbers of gene trees supporting two comparable duplication events at _s_<br>
**Note, to obtain Total MEscore value for the species tree, one should sum all values in ME score column.**

## Examples:

`python gdilp.py tests/treefam.txt ` <br>
Computes GME model (called also unrestricted, FHS), that is, NP-hard Problem 1.The output is written to out_treefam_0_-1.csv.<br>

`python gdilp.py --intervallimit 1 tests/treefam.txt`<br>
Computes 0-GME model (LCA). The output is written to out_treefam_1_-1.csv.<br>

`python gdilp.py --max_converted_spec 0 tests/treefam.txt`<br>
Computes &infin;,0-GME model (RME). The output is written to out_treefam_0_0.csv.<br>

`python3 gdilp.py --intervallimit 7 --max_converted_spec 3066 --time_limit 25000 tests/treefam.txt` <br>
Computes 7,3066-GME model. The output is written to out_treefam_7_3066.csv. <br>
If the computation not ends in 25000 seconds, then the output consist the current best solution. <br>


## Publication:

Novel genomic duplication models through integer linear programming
 
Jarosław Paszek, Oliver Eulenstein, Paweł Górecki

BCB '21: Proceedings of the 12th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics<br> 
August 2021; 
Article No.: 16 Pages 1–11 https://doi.org/10.1145/3459930.3469549



