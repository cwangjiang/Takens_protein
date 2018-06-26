# Takens_protein
This documentation is for applying delay embedding on proteins.

# 1 Trp-cage (2JOF): 

## 0_simulation

There are two trajectories, the smaller one contains only alpha-C atoms, the larger one contains all protein atoms, we are using the smaller alpha-C one, extract it.

The long simulation is splited into 10 .dcd files, load them all together into VMD:

```bash
vmd 2JOF-0-c-alpha.mae 2JOF-0-c-alpha-000.dcd 2JOF-0-c-alpha-001.dcd 2JOF-0-c-alpha-002.dcd 2JOF-0-c-alpha-003.dcd 2JOF-0-c-alpha-004.dcd 2JOF-0-c-alpha-005.dcd 2JOF-0-c-alpha-006.dcd 2JOF-0-c-alpha-007.dcd 2JOF-0-c-alpha-008.dcd 2JOF-0-c-alpha-009.dcd 2JOF-0-c-alpha-010.dcd
```

Then use VMD to transfer the trajector into .pdb corrdinates `trajectory_all_stride1.pdb`.
Stride1 means there is no subsample, the interval between two frame is 0.2ns, and there are 1,000,000 points


## 1_dMaps

- Move `trajectory_all_stride1.pdb` here, and transfer .pdb to .gro file:
```bash
editconf -f trajectory_all_stride1.pdb -o traj_pbc.gro
```
This will generate `traj_pbc.gro`.

- Use `clean_data.cpp` to clean .gro file, making it more structured:
```bash
g++ clean_data.cpp -o cleandata.out

./cleandata.out
```
this will generate `traj_Calpha.dat`.

- Load this structured `traj_Calpha.dat` into matlab, save it as an -ascii format .mat file: `traj_Calpha.mat`, so that it can be loaded in C++.
```bash
>> load traj_Calpha.dat

>> save('traj_Calpha.mat','traj_Calpha','-ascii')
```

- use `subtraj.m` in matlab to subsample 100,000 (skip size 10), or 10,000 (skip size 100) points, make it smaller to handle, we will use 100,000 points, because 1 million points are too much, 100k is OK, but still slow to run multiple times for finding parameters and debuging, so we will use the smaller set with 10,000 points to find hyper parameters, and then use such parameters to apply to 100,000 trajectory. 
```bash
>> subtraj.m
```
This will generate big sampling `traj_Calpha_skip10.mat` and smaller one `traj_Calpha_skip100.mat`.

- Then use `main.cpp` to compute pairwise distances and select out pivots for pivot-diffusion maps. First use small set with 10,000 to find good rcut, in `main.cpp`, set N = 10000, numcut = 100, trial rcut = 0.41, traj.load("traj_Calpha_skip100.mat",raw_ascii), then execute:
```bash
g++ -std=c++0x  main.cpp functions.cpp -o main.out -O2 -larmadillo -llapack -lblas

./main.out
```
This will generate `Distance.mat`, `pivot.mat`, `pivotindex.mat`. 

- Use `loadingfile.m` in matlab to load these matrixs, and find out the number of pivots.
	- `Distance.mat` is the m by N distance matrix, where N is total point number 10,000, m is the number of pivots. 
	- `pivot.mat` is a m by 3 matrix also provides information of the pivots: the first column is index 1 to m, the second column is the index of the point form all N size that this pivot corresponds to, the third column is the domain size of this pivot, which is bounded by the number cut off numcut, which should not be too large, otherwise such pivot will include too many point and desory the local fine structures, empirically this hyper parameters is set to be N\*0.01~N\*0.1. 
	- `pivotindex.mat` is list out all N points and the pivot ID it belongs to.
```bash
>> loadingfile
```

We find that using rcut = 0.41 give us 400 pivots, usually around 500~1000 pivots is pretty enough to cover the whole manifold, if we apply these parameters on the larger 100,000 points set, it should also give us 400 points, because these 400 points already cover the manifold, but in practice, we find it will give us more pivots on larger data sets, this should be due to larger data sets samples more unexplored regions. Then we change N = 100000, Nct = 1000, rcut = 0.41, traj.load("traj_Calpha_skip10.mat",raw_ascii), and excute the above codes, we get 622 pivots from the 100,000 points set.

- `loadingfile.m` will generate the 622 by 622 distance matrix of only pivots, then use `dMap.m` in matlab to conduct diffusion maps on 622 pivot points matrix. In the `dMap.m` code, we need set N = 622, then tune parameters eps and <img src="https://latex.codecogs.com/gif.latex?\alpha">, where <img src="https://latex.codecogs.com/gif.latex?\alpha"> is the factor to rescale all pairwise distances, to balance the density, usually <img src="https://latex.codecogs.com/gif.latex?0.1<\alpha<1.0">, we can try different one, and see if it can shrink the sparse and expaned the condensed region to reveal more fine structures, we can also compute local density, to make sure good <img src="https://latex.codecogs.com/gif.latex?\alpha"> make the difference between maximum and minimum local density reasonable. eps is the Gaussian kernal bandwidth, it should be larger than rcut. We set esp = 1, and 
<img src="https://latex.codecogs.com/gif.latex?\alpha=0.15">, then execute:
```bash
>> dMap
```
This will generate `dMap.mat` which store the eigenvectors and eigenvalues for 622 pivots. We find that top two non-trivial eigenvectors are <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">.

- Then use `nystrom.m` in matlab to insert the rest 100,000-622 points back in to the diffusion maps, set N = 100000, esp = 1, <img src="https://latex.codecogs.com/gif.latex?\alpha=0.15">, so that they are consistent with the dMaps:
```bash
>> nystrom
```
This generate `X.mat`, which is an N by 9 matrix, and the embedding of all 100,000 points in to the dMap space, but we only care about the first two colums, since they correspond to <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">.


- Then compute head to tail distance using `compute_h2t.m`, here we don't use the subsamplings, we use the full trajectory with 1,000,000 points, so that the resolution of the h2t(t) is high enough.
```bash
>> compute_h2t
```
this will generate `h2t.mat`.

- Use `compute_RMSD.m` to compute Rg, RMSD, RMSD_helix, the native state is stored as `trpnativeCalpha.mat`
```bash
>> compute_RMSD
```
This will generate `RMSD.mat`, `RMSD_helix.mat`, and `Rg.mat`, Rg is a N by 4 matrix, the first column is Radius of gyration, the 2,3,4 columns are 1st, 2nd, 3rd component of gyration tensor. 

- Then we can use `plot_Rg.m` to plot dmaps, and each point is colored as Rg/RMSD/RMSD_helix.
```bash
>> plot_Rg
```

- Then use `FES.m` to compute FES in diffusion map space <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">.
```bash
>> FES
```

- Then use `FES_new.m` to compute FE associate with each point, and can be used for computing free energy correlation with the reconstructed one. 
```bash
>> FES_new
```
this will generate `OrigionFE.mat`.

- Use `FES_RMSD.m` to compute FE and plot FE in the conventional space.
```bash
>> FES_RMSD
```

## 2_MI
Move `h2t.mat` from 1_dMaps here, and use `MI.m` to compute mutual information. Good delay time is where MI decays to 1/e.
```bash
>> MI
```
We select the bin to be 5 for running MI, and set the maximum tau to be 100, which equal to 20ns, because the interval between two adjacent points is 0.2ns. We find that 0.6ns is a good delay time. 

## 3_FNN
Move `h2t.mat` here, and use `FNN.m` or `FNN_fast.m` to compute E1 function, which give good delay dimension D, `FNN_fast.m` is optimized from `FNN.m`, may runs faster.
```bash
>> FNN_fast(h2t)
```
E1 saturate at about 20 dimension, we can use any dimension greater than 20, and we use 50 dimension.

## 4_RCT
Move `h2t.mat` here, and use delay time tau and delay dimension D to generate a reconstucted euclidean space, `EBD.mat`
```bash
>> RCT
```
This will generate a 100,000 by 50 matrix, each row is a reconstructed point, there are 100,000 points, the full length of `h2t.mat` is 1,000,000, and we sub-sample 100,000 from it, with a skip size of 10. 

## 5_RCT_dMaps
Similar to 1\_dMaps, we apply pivot diffusion maps to the delayed points `EBD.mat`, distance between two point is not RMSD between molecular configurations, but regular Eucledian distance. The brief process and parameters are provided below, details can be referred to 1_dMaps.

- Use `main.cpp` to compute pairwise distances, setting numcut = 10000, and rcut = 8.5, traj.load("EBD.mat",raw_ascii), it will output 515 pivots. Similar to 1_dMaps, since 100,000 points are too much and runs slow, we can test different rcut on a smaller sampling of 10,000 points `EBDsub.mat` first, then apply the rcut to the big trajectory `EBD.mat`.

- Then using `loadingfile.m` to load the pairwise distance into matlab. 

- Use `dMap.m` to conduct diffusion maps on the pivots, where eps = 2.1 and alpha = 0.15.

- Use `nystrom.m` to insert all points, and generate `X.mat`.

- Use `plot_Rg.m` to plot reconstructed dmaps, and each point is colored as Rg/RMSD/RMSD_helix.

- Use `FES.m` to plot FES on diffusion map space <img src="https://latex.codecogs.com/gif.latex?\psi_2^*,\psi_4^*">. Since <img src="https://latex.codecogs.com/gif.latex?\psi_2^*,\psi_3^*"> are functional correlated, so we ignore <img src="https://latex.codecogs.com/gif.latex?\psi_3^*">.

- Use `FES_new.m` to compute FE for each point, and generate `delayFE.mat`. Move `OriginFE.mat` here, and use `Compute_Correlation.m` to compute free energy correlation.
```bash
>> Compute_Correlation
```

## 6_detJ
Move `X.mat` from 1\_dMaps and rename it as `X_origin.mat`, move `X.mat` from 5_RCT_dMaps and rename it as `X_delay.mat`. Use `/bandwidthscan/bandwidthscan.m` to find correct Gaussian kernel bandwidth for computing detJ. We find that the measurement of the error decrease to 1/e at around 0.3. Then use `plotDetj.m` to compute detJ for the first 10000 points:
```bash
>> plotDetj
```
In the plotDetj.m, change `meshless_jacobian(X_origin,X_delay,0)` into `meshless_jacobian(X_delay,X_origin,0)` to compute the detJ for the backward mapping.


# 2 Villin: 

2 Villin is similar to 1 Trp-cage, only parameters are provided below:

- Simulation is 120,000ns, the interval between two samplings is 0.2ns, totally 600,000 points, we subsample 60,000 to perform analysis.

- In the pivot diffusion maps of all atom configuration, numcut = 1000, rcut = 0.8, which gives 448 pivots, in diffusion maps, eps = 1.0, alpha = 0.7, then extract CV <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">.

- In the delay embedding, delay time is 0.2*2 = 0.4ns, delay dimension D is 50. 

- In the diffusion maps on the reconstructed data, numcut = 6000, rcut = 8.0, this gives 1358 pivots, and eps = 9, alpha = 0.5, and extract two top CV <img src="https://latex.codecogs.com/gif.latex?\psi_2^*,\psi_3^*">.

- The original well depth is 3.76, reconstructed well depth is 3.90, and the FE correlation is 0.46. 

- When computing Jacibian determinant, bandwidth is selected as 0.3.



# 3 BBA: 

3 BBA is also similar to 1 Trp-cage, only parameters are provided below:

- Simulation is 200,000ns, the interval between two samplings is 0.2ns, totally 1,000,000 points, we subsample 100,000 to perform analysis.

- In the pivot diffusion maps of all atom configuration, numcut = 1000, rcut = 1.0, which gives 596 pivots, in diffusion maps, eps = 1.1, alpha = 0.15, then extract CV <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">.

- In the delay embedding, delay time is 0.2*4 = 0.8ns, delay dimension D is 50. 

- In the diffusion maps on the reconstructed data, numcut = 10000, rcut = 8.0, this gives 840 pivots, and eps = 10, alpha = 0.5, and extract two top CV <img src="https://latex.codecogs.com/gif.latex?\psi_2^*,\psi_3^*">.

- The original well depth is 2.60, reconstructed well depth is 3.77, and the FE correlation is 0.56. 

- When computing Jacibian determinant, bandwidth is selected as 0.5.




# 4 Trp-cage (1L2Y): 
This part deals with 37 individual simulation, and will be slightly different from 1 Trp-cage:

## 0_Simulation
Full simulation trajectories of four types of simulation are provided: 1. wild type Trp-cage at different temperatures, 2. designed mutation (in Nature paper https://www.nature.com/articles/nsb798), 3. single Alanine scanning mutations, 4. tetrad Alanine scanning mutations. Mutated configuration is generated throught PEP-FOLD3: http://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3/, for problematic configurations, for example, there is missing residue or missing atoms, use PDBFixer to fix it: https://github.com/pandegroup/pdbfixer. Then run the simulation using openMM, the input and output files for each simulation are provided. /Cluster contains all full simulations from cluster.

Each simulation is 1 micro second, containing 100,000 points, with interval to be 0.01ns. For each simulation from openMM, the trajectory is in `output.pdb`, we first need to generate a coordinate reference from the .pdb file using editconf, and then transfer the long `output.pdb` into .gro file, which only contain 20 alpha-carbon atoms, then clean the .gro file and subsample 10,000 points. All these are include in the `CreatCoord.sh` script.
```bash
./CreatCoord.sh
```
This will generate `coordinates_sub_XXX.dat` file the subsampling of only 10,000 configurations, and each configuration contain only alpha-carbon atoms.

## 1_dMaps

#### 1_All_subtraj
Move all coordinates_sub_XXX.dat here, use `combinetraj.m` to combine all trajectories into one ensemble in matlab.
```bash
>>Combinetraj
```
This will generate `traj.mat`, which contain 370,000 points.
We also generate a ensemble even smaller trajectory, each simulation just contatin 1000 point, the ensemble is `subtraj.mat`, which contain 37,000 points.

#### 2_pdmap
Move `subtraj.mat` here to conduct pivot diffusion maps on 37,000 points. Similar as before, numcut = 1000, rcut = 0.38, which gives 478 pivots, in diffusion maps, eps = 0.8, alpha = 1.0, then extract CV <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">. This will generate a `dMap.mat`, which is the diffusion map of the small ensemble. 

#### 3_nystrom
Move `dMap.mat` `traj.mat` here, we use  `main.cpp` to compute distances between all 370,000 points and the 478 pivots, and insert 370,000 points back into the <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3"> space. This will generate `X.mat`.

In /FES, use `RUNFES.m` to compute and plot all FES for each of 37 systems. `RUNFES_png.m` is to plot and save a .png file. `RUNFES_new.m` is to compute free energy associate with each point for each system, and generate `fes.mat`.

In /FES/Differences, use use `RUNFES.m` to compute and plot all FES difference for each of 37 systems relative to #5, the reference FES is stored as Href.mat. `RUNFES_png.m` is to plot and save a .png file. `RUNFES_new.m` is to compute free energy differences associate with each point for each system, and generate `fesD.mat`.

#### 4_conventional

Similar to 3_nystrom, but compute and plot all FES in the conventional space. 

## 2_RCT
To conduct delay embedding and diffusion maps on the reconstructed points. 

#### 1_MI
To compute delay time using the fifth simulation, which is the wild type at 380K, because the simulation are indexing from 0, it show `h2t_4.mat`, not 5, which is a little bit weird. Delay time is determined to be 0.02ns.

#### 2_FNN
FNN determines the delay dimension to be 20.
```bash
>>FNN_fast(ht2_4)
```

#### 3_RCT
Construct delayed points using the above delay time and dimension for each simulation individually, and then combine all reconstructed points into a bigger ensemble `traj10k.mat`, where each simulation contribute 10k delayed points, and `traj1000.mat`, where each simulation contributes 1000 points. 

#### 4_pdMap
Apply `main.cpp` on the smaller ensemble traj1000.mat, with numcut = 2000, rcut = 2.5, this identify 970 pivots, then apply diffusion maps on these 970 pivots with eps = 4.0 and alpha = 0.25. then extract CV <img src="https://latex.codecogs.com/gif.latex?\psi_2^*,\psi_4^*">.  This generates dMap.mat.

#### 5_nystrom
Compute distance between all 370,000 reconstructed point to the 970 pivots, and use nystrom to insert all 370,000 points back into the manifold, which generate X.mat. 

In /FES, use `RUNFES.m` to compute and plot all FES for each of 37 systems. `RUNFES_png.m` is to plot and save a .png file. `RUNFES_new.m` is to compute free energy associate with each point for each system, and generate `fesR.mat`.

In /FES/Differences, use use `RUNFES.m` to compute and plot all FES difference for each of 37 systems relative to #5, the reference FES is stored as Href.mat. `RUNFES_png.m` is to plot and save a .png file. `RUNFES_new.m` is to compute free energy differences associate with each point for each system, and generate `fesRD.mat`.

Move `fes.mat` here and use `Compute_corr.m` to compute the FE correlation between the original FES and the reconstructed FES for each of the simulation. And compute the well depth for the original and the reconstructed one. 

## 3_detJ

Move `X.mat` from 1_dMaps/3_nystrom/X.mat here and rename it `X_Original.mat`. Move `X.mat` from 2_RCT/5_nystrom/X.mat and rename it `X_delay.mat`. Use `plotDetj.m` to compute Jacobian determinant.














