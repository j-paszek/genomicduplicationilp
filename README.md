# Usage: 
<pre>
`gdilp.py [-h] [--intervallimit [INTERVAL_LIMIT]]
         [--max_converted_spec [MAX_CONVERTED_SPEC]]
         [--time_limit [TIME_LIMIT]]
         input_file`
 </pre> 

## positional arguments:

`input_file`  &ensp;Default input is a file in RME format. The first line consist of the number of gene trees n. <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; The next n lines represent gene trees. Then, the following line consist of the species tree.<br>

## optional arguments:

`-h, --help`  &ensp; show help message and exit<br>

`--intervallimit [INTERVAL_LIMIT]` Set the maximal length of intervals.<br> 
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;  Default k=0, and no intervals are shortened. It is equivalent to $\lambda = \infty$ used in the publication.<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; If equals k>0, then all intervals with length y>k are shortened to the size k. <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Note &lambda; $\lambda$ from publication (edge interval length) equals $k-1$ (as parameter k denote node length of an interval)<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; In example for a gene tree node $g$ the interval equals $\langle M(g),b \rangle$, where $M(g)$ is a lca-mapping of $g$ in the species tree,<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;  and $b$ is a node in species tree. If $|\langle M(g),b \rangle| = y > k$ then we assign new interval $\langle M(g),c \rangle$ such that <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; $|\langle M(g),b \rangle| = k$ and $c \in \langle M(g),b \rangle$;  nodes are removed from 'b' side. <br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; In other words we limit the distance between the duplication mapping and lca-mapping not to exceed k.<br>

`--max_converted_spec [MAX_CONVERTED_SPEC]` Set the maximal number of speciations that can be converted into duplications ($\sigma$)<br>
&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Default parameter $\sigma = \infty$. To manually set $\sigma = \infty$, write -1. <br>

`--time_limit [TIME_LIMIT]` Set the maximal computation time in seconds.<br>



# Installation:

**Requirements**

(A) Installation gurobi solver.<br> 
- _Local machines_ (A) be sufficient to run `gdilp` (install also Biopython, Pandas if not present)<br>
- _Remote servers_ Please try (B)<br>
(B) Create python environment as described below.<br>

Default modules: &emsp; &emsp;&emsp;&ensp; timeit,  argparse (Python >= 3.2), io<br>
Required modules: &emsp; &emsp; &ensp;Bio, gurobipy, pandas<br>
Optional:&emsp; &emsp;&emsp;&emsp; &emsp;&emsp; &emsp;&ensp;pytest<br>

Info:<br>
- See `requiremets.txt` for other requirements info (eg. gurobipy version 9.1.1 was used).<br>
- The installation of gurobipy requires sudo, therefore, I suggest to create a python environment as described below.<br>

**Download sources**<br>

Download zipped folder of this repository. Unzip to project source folder (like `genomicduplicationilp-master/`).<br>

**Create python environment**<br>

(1) `python3 -m venv gurobi-env` &emsp;&emsp; &emsp;&emsp;&emsp;&emsp; Create virtual environment at project source folder.<br>
(2) `source ./gurobi-env/bin/activate` &emsp;&emsp; &emsp;Activate the environment.<br>
gurobipy installation:<br>
(3) `cd /opt/gurobi901/linux64/` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp; Go to the gurobi folder.<br>
(4) `python3 setup.py build -b PATH install` Gurobi installation to a selected PATH. Useful on remote servers. <br>
&emsp;&emsp;Running only this command not in virtual environment is not enough, because it instals gurobipy into the python folder <br>
&emsp;&emsp;causing another permission denied (not in /opt/gurobi... but in python folder instead)<br>
(5) `cd genomicduplicationilp-master/` &emsp;&emsp;&emsp;&ensp; Go to the project source folder.<br>
(6) `pip3 install biopython ` &emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&emsp; &emsp;&emsp;&ensp; Biopython installation.<br>
(7) `pip3 install pandas` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp;&ensp; Pandas installation.<br>
Now the environment is ready to use. Use commands (i)-(iv). Than multiple runs like (v). See below.<br>
(8) `deactivate` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Close virtual environment.<br>

**Run from environment**

( i ) `source ./gurobi-env/bin/activate` &emsp;&emsp; &emsp;&emsp;&emsp;&emsp;Activate the environment.<br>
(ii )  `export GUROBI_HOME="/opt/gurobi901/linux64"` Now 3 magic configuration lines.<br>
(iii)  `export PATH="${PATH}:${GUROBI_HOME}/bin"`<br>
(iv) `export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"`<br>
(v ) `python3 gdilp.py tests/guigo.txt` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Run the ILP.<br>
(vi)  `deactivate`&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; &emsp;&emsp;  Close virtual environment.<br>

Note:<br>
Normally, we only need to activate the environment. And than run as many cumputations (v) as we desire.<br>
However, if we log out from the server and then log in again, then we have to run 3 magic lines again.<br>

## Project files

**Simulated data**<br>
The original SimPhy results are in `wgd-simulations` folder.<br>
`sim_utils.py` was used to extract data to input files stored as `tests\sim[i].txt`<br>
Note `[i]`file corresponds to $S_i$ tree (see Figure 1)<br>

**Input data**<br>
Input data is stored in `tests\` folder. In our study we use TreeFam dataset - file `treefam.txt`.<br>
The folder contains also sample tests and `test_gd_ilp.py` to perform unittest.<br>

**Output data**<br>
The experimental section describes the evaluation of the results stored at `results\out_name_k_`$\sigma$`.csv`<br>
Where:<br>
- name - sim1, ..., sim5, treefam - is the name of the corresponding input dataset<br>
- k - the value of INTERVAL_LIMIT parameter used<br>
- $\sigma$ -the value of MAX_CONVERTED_SPEC parameter used<br>

## Examples:

`python gdilp.py tests/treefam.txt ` <br>
Computes GME model (called also unrestricted, FHS), that is, NP-hard Problem 1.The output is written to out_treefam_0_-1.csv.<br>

`python gdilp.py --intervallimit 1 tests/treefam.txt`<br>
Computes 0-GME model (LCA). The output is written to out_treefam_1_-1.csv.<br>

`python gdilp.py --max_converted_spec 0 tests/treefam.txt`<br>
Computes $\infty$,0-GME model (RME). The output is written to out_treefam_0_0.csv.<br>

`python3 gdilp.py --intervallimit 7 --max_converted_spec 3066 --time_limit 25000 tests/treefam.txt` <br>
Computes 7,3066-GME model. The output is written to out_treefam_7_3066.csv. <br>
If the computation not ends in 25000, then the output consist the current best solution. <br>



