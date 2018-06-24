# Takens_protein
The documentation for applying delay embedding on proteins

# 1 Trp-cage (2JOF): 

## 0_simulation

There are two trajectories, the smaller one contain only alpha-C atoms, the larger one contain all protein atoms, we are using the smaller alpha-C one, extract it.

The long simulation is splited into 10 .dcd files, load them all together into VMD:

```bash
vmd 2JOF-0-c-alpha.mae 2JOF-0-c-alpha-000.dcd 2JOF-0-c-alpha-001.dcd 2JOF-0-c-alpha-002.dcd 2JOF-0-c-alpha-003.dcd 2JOF-0-c-alpha-004.dcd 2JOF-0-c-alpha-005.dcd 2JOF-0-c-alpha-006.dcd 2JOF-0-c-alpha-007.dcd 2JOF-0-c-alpha-008.dcd 2JOF-0-c-alpha-009.dcd 2JOF-0-c-alpha-010.dcd
```

Then use VMD to transfer the trajector into .pdb corrdinates. trajectory_all_stride1.pdb
Stride1 means there is no subsample, the interval between two frame is 0.2ns, and there are 1,000,000 points


## 1_dMaps

- To transfer .pdb to .gro file:
```bash
editconf -f trajectory_all_stride1.pdb -o traj_pbc.gro
```
This will generate traj_pbc.gro

- Use `clean_data.cpp` to clean .gro file, making it more structured:
```bash
g++ clean_data.cpp -o cleandata.out

./cleandata.out
```
this will generate traj_Calpha.dat

- Load this structured traj_Calpha.dat into matlab, save it as an -ascii format .mat file: traj_Calpha.mat, so that it can be loaded in C++.
```bash
>> load traj_Calpha.dat

>> save('traj_Calpha_test.mat','traj_Calpha','-ascii')
```

- use `subtraj.m` in matlab to subsample 100,000 (skip size 10), or 10,000 (skip size 100) points, make it smaller to handle, we will use 100,000 points. because 1 million is too much, 100k is OK, but still slow to run multiple times for finding parameters and debuging, so we will use the smaller set with 10,000 points to find hyper parameters, and then use such parameters to apply to 100,000 trajectory. 
```bash
>> subtraj.m
```

- Then use `main.cpp` to compute pairwise distances and select out pivots for pivot-diffusion maps. First use small set with 10,000 to find good rcut, in main.cpp, set N = 10000, Ncut = 100, trial rcut = 0.41, traj.load("traj_Calpha_skip100.mat",raw_ascii), then execute:
```bash
g++ -std=c++0x  main.cpp functions.cpp -o main.out -O2 -larmadillo -llapack -lblas

./main.out
```
This will generate Distance.mat, pivot.mat, pivotindex.mat. 

- Use `loadingfile.m` in matlab to load these matrix, and find out the number of pivots.
	- Distance.mat is the m by N distance matrix, where N is total point number 10,000, m is the number of pivots. 
	- pivot.mat is a m by 3 matrix also provides information of the pivots: the first column is index 1 to m, the second column is the index of the point form all N size that this pivot corresponds to, the third column is the domain size of this pivot, which is bounded by the number cut off Ncut, which should not be too large, otherwise such pivot will include too many point and desory the local find structures, empirically this hyper parameters is set to be N\*0.01~N\*0.1. 
	- pivotindex.mat is list out all N points and the pivot ID the belong to.
```bash
>> loadingfile
```

We find that using rcut = 0.41 give us 400 pivots, usually around 500~1000 pivots is pretty enough to cover the whole manifold, if we apply these parameters on the larger 100,000 points set, it should also give us 400 points, because these 400 points already cover the manifold, but in practice, we find it will gives us more pivots on larger data sets, this should be due to larger data sets samples more unexplored regions. Then we change N = 100000, Nct = 1000, rcut = 0.41, traj.load("traj_Calpha_skip10.mat",raw_ascii), and excute the above codes, we get 622 pivots from the 100,000 points set.

- loadingfile.m will generate the 622 by 622 distance matrix of only pivots, then use dMap.m in matlab to conduct diffusion maps on 622 pivot points matrix. In the dMap.m code, we need set N = 622, then tune parameters eps and <img src="https://latex.codecogs.com/gif.latex?\alpha">, where <img src="https://latex.codecogs.com/gif.latex?\alpha"> is the factor to rescale all pairwise distances, to balance the density, usually <img src="https://latex.codecogs.com/gif.latex?0.1<\alpha<1.0">, we can try different one, and see if it can shrink the sparse and expaned the condensed region to reveal more fine structures, we can also compute local density, to make sure good <img src="https://latex.codecogs.com/gif.latex?\alpha"> make the difference between maximum and minimum local density reasonable. eps is the Gaussian kernal bandwidth, it should be larger than rcut. We set esp = 1, and 
<img src="https://latex.codecogs.com/gif.latex?\alpha=0.15">, then execute:
```bash
>> dMap
```
This will generate dMap.mat which store the eigenvectors and eigenvalues for 622 pivots. We find that top two non-trivial eigenvectors are <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">.

- Then use `nystrom.m` in matlab to insert the rest 100,000-622 points back in to the diffusion maps, set N = 100000, esp = 1, <img src="https://latex.codecogs.com/gif.latex?\alpha=0.15">, so that they are consistent with the dMaps:
```bash
>> nystrom
```
This generate X.mat, which is an N by 9 matrix, and the embedding of all 100,000 points in to the dMap space, but we only care about the first two colums, since they correspond to <img src="https://latex.codecogs.com/gif.latex?\psi_2,\psi_3">.


- Then compute head to tail distance using compute_h2t.m, here we don't use the subsamplings, we use the full trajectory, so that the resolution of the h2t(t) is high enough.
```bash
>> compute_h2t
```
this will generate h2t.mat.

- Use `compute_RMSD.m` to compute Rg, RMSD, RMSD_helix, the native state is stored as trpnativeCalpha.mat
```bash
>> compute_RMSD
```
This will generate RMSD.mat, RMSD_helix.mat, and Rg.mat, Rg is a N by 4 matrix, the first column is Radius of gyration, the 2,3,4 columns are 1st, 2nd, 3rd component of gyration tensor. 

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
>> FES
```
this will generate OrigionFE.mat

- Use `FES_RMSD.m` to compute FE and plot FE in the conventional space.
```bash
>> FES_RMSD
```

## 2_MI
Move h2t.mat from 1_dMaps here, and use `MI.m` to compute mutual information. good delay time is where MI decays to 1/e.
```bash
>> MI
```
We select the bin to be 5 for running MI, and set the maximum tau to be 100, which equal to 20ns, because the interval between two adjacent points is 0.2ns. We find that 0.6ns is a good delay time. 

## 3_FNN
Move h2t.mat here, and use FNN.m or `FNN_fast.m` to compute E1 function, which give good delay dimension D, FNN_fast.m is optimized from FNN.m, may runs faster.
```bash
>> FNN_fast(h2t)
```
E1 saturate at about 20 dimension, we can use any dimension greater than 20, and we use 50 dimension.

## 4_RCT
Move h2t here, and use delay time tau and delay dimension D to generate a reconstucted euclidean space, EBD.mat
```bash
>> RCT
```
This will generate a 100,000 by 50 matrix, each row is a reconstructed point, there are 100,000 points, the full length of h2t.mat is 1,000,000, and we sub-sample 100,000 from it, with a skip size of 10. 

## 5_RCT_dMaps
Similar to 1_dMaps, we use pivot diffusion maps to the delayed points EBD.mat, distance between two point is not RMSD between molecular configuration, but regular Eucledian distance. The brief process and parameters are provided below, details can be referred to 1_dMaps.

- Use `main.cpp` to compute pairwise distances, setting Ncut = 10000, and rcut = 8.5, traj.load("EBD.mat",raw_ascii), it will output 515 pivots. Similar to 1_dMaps, since 100,000 points are too much and runs slow, we can test different rcut on a smaller sampling of 10,000 points EBDsub.mat firt, then apply the rcut to the big trajectory EBD.mat.

- Then using `loadingfile.m` to load the pairwise distance into matlab. 

- Use `dMap.m` to conduct diffusion maps on the pivots, where eps = 2.1 and alpha = 0.15.

- Use `nystrom.m` to insert all points, and generate `X.mat`.

- Use `plot_Rg.m` to plot reconstructed dmaps, and each point is colored as Rg/RMSD/RMSD_helix.

- Use `FES.m` to plot FES on diffusion map space <img src="https://latex.codecogs.com/gif.latex?\psi_2^*,\psi_4^*">. Since <img src="https://latex.codecogs.com/gif.latex?\psi_2^*,\psi_3^*"> are functional correlated, so we ignore <img src="https://latex.codecogs.com/gif.latex?\psi_3^*">..

- Use `FES_new.m` to compute FE for each point, and generate `delayFE.mat`. Move OriginFE.mat here, an use `Compute_Correlation.m` to compute free energy correlation.
```bash
>> Compute_Correlation
```

- Use `plot_Correlation.m to compute the FE correlation between original and reconstructed FES. 

## 6_detJ
Move X.mat from 1_dMaps and rename it as X_Original.mat, move X.mat from 5_RCT_dMaps and rename it as X_delay.mat. Use `/bandwidthscan/bandwidthscan.m` to find correct Gaussian kernel bandwidth for compute detJ. We find that the measurement of the error decrease to 1/e at around 0.3. Then use `plotDetj.m` to compute detJ for the first 10000 points:
```bash
>> plotDetj
```
In the plotDetj.m, change `meshless_jacobian(X_origion,X_delay,0)` into `meshless_jacobian(X_delay,X_origin,0)` to compute the detJ for the backward mapping.


# 2 Villin: 

2 Villin is similar to 1 Trp-cage



# 3 BBA: 

3 BBA is similar to 1 Trp-cage




# 4 Trp-cage (1L2Y): 
