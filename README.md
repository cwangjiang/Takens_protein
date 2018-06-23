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

to transfer .pdb to .gro file:
```bash
editconf -f trajectory_all_stride1.pdb -o traj_pbc.gro
```
This will generate traj_pbc.gro

Use clean_data.cpp to clean .gro file, making it more structured:
```bash
g++ clean_data.cpp -o cleandata.out
./cleandata.out
```
this will generate traj_Calpha.dat

Load this structured traj_Calpha.dat into matlab, save it as a .mat file: traj_Calpha.mat

use subtraj.m in matlab to subsample 100,000 (skip10), or 10,000 (skip100) points, make it smaller to handle, we will use 100,000 points. because 1 million is too much, 100k is OK, but still slow to run multiple time for finding parameters, so skip 100, with 10,000 points small sampling is to find parameters, and then use such parameters to apply to 100,000 trajectory. 

Then use main.cpp to conduct pivot-diffusion maps, here is just to compute the pairwise distances:
```bash
g++ -std=c++0x  main.cpp functions.cpp -o main.out -O2 -larmadillo -llapack -lblas
./main.out
```
This will generate Distance.mat/pivot.mat/pivotindex.mat Distance is the N by m distance matrix, N is 100,000, m is the number of pivots

Use loadingfile.m in matlab to load these matrix, and find out the number of pivots.

Use dMap.m in matlab to conduct diffusion maps on the m by m pivot points matrix. 

Then use nystrom.m in matlab to insert the rest N-m points back in to the diffusion maps, and get the result matrix: X.mat

Notice: In main.cpp, we loop over all 100,000 points to find out about 500-1000 pivots, which depend on the number cutoff and radius cutoff of each domain, the criterial is: Ncut is about 1-10% of the totall N, and test different rcut so that it generate about 500 pivot for you, but the test is slow on 100,000 points, so we can do this test on a smaller smaplling, once we know the good rcut, we can apply it to the 100,000 points. In dMaps, we also need a \epsilon in Gaussian, and \alpha. the criterial is: \epsilong should be slightly larger than rcut. Then you compute the local density of each point using Density.m, compute the ration of the max_density/min_density, this may be a large number, for example say R = 10^50, \alpha should rescale all distance and make R smaller: R^(\alpha) = 10^5, so \alpha = 0.1.

Then compute head to tail distance using compute_angle.m 
this will generate distance120.mat

Use compute_RMSD.m to compute Rg,RMSD, RMSD_helix, the native state is stored  as trpnativeCalpha.m
This will generate RMSD.mat

Then we can use plot_Rg.m to plot dmaps, and each point is colored as Rg/RMSD.

Then use FES.m to compute FES in the diffusion map space.

Then use FES_new.m to compute FE associate with each point, this will generate OrigionFE.mat

Use FES_RMSD.m to compute FE and plot FE in the conventional space.

## 2_MI
use MI.m to compute mutual information. good delay time is where MI decays to 1/e.

## 3_FNN
Use FNN.m or FNN_fast.m to compute E1 function, FNN_fast is a optimized verstion, may runs faster. 

## 4_RCT
Use delay time and delay dimension to generate a reconstucted euclidean space, EBD.mat

## 5_RCT_dMaps
Similar to dMaps, we use pivot diffusion maps to the delayed points EBD.mat, distance between point is not RMSD between molecular configuration, but regular Eucledian distance.

Use 
  
main.cpp to compute pairwise distances, and generate about 500 pivots, if 100,000 is too much and runs slow, we can try this on a smaller sampling of 10,000 points. Determine the Ncut first, and try different rcut. 

Then using loadingfile.m to load the pairwise distance into matlab. 

Use dMap.m to conduct diffusion maps on the pivots.

Use nystrom.m to insert all points, and generate X.mat.

Use FES.m to plot FES on diffusion map space

Use FES_new.m to compute FE for each point, and generate delayFE.mat

Use plot_Correlation.m to compute the FE correlation between original and reconstructed FES. 

## 6_detJ
Use meshless_jacobian.m to compute detJ between original and reconstructed FES. 




# 2 Villin: 




# 3 BBA: 




# 4 Trp-cage (1L2Y): 
