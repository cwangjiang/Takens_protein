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

## 2_MI

## 3_FNN

## 4_RCT

## 5_RCT_dMaps

## 6_detJ



# 2 Villin: 




# 3 BBA: 




# 4 Trp-cage (1L2Y): 
