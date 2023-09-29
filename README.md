# Approximate Component Mode Synthesis (ACMS) method for heterogeneous Helmholtz equations

This repository contains the Matlab implementation of the extension of the ACMS method to heterogeneous Helmholtz equations as introduced in the paper

**An extension of the approximate component mode synthesis method to the heterogeneous Helmholtz equation** <br>
_Elena Giammatteo, Alexander Heinlein, Matthias Schlottbom_ <br>
[https://arxiv.org/abs/2303.06671](https://arxiv.org/abs/2303.06671)

The Matlab scripts for the implementation of the method can be found in the folder `functions`, and the scripts in the main directory correspond to the examples tested in the paper:
+ `run_boundary_source.m` - Classical Helmholtz example with boundary source; see _Section 5.3_.
+ `run_localized_interior_source.m` - Classical Helmholtz example with boundary source, see _Section 5.2_.
+ `run_periodic.m` - Periodic structure example with planewave solution; see _Section 5.4_.
+ `run_planewave.m` - Classical Helmholtz example with planewave solution; see _Section 5.1_.

If you use this code, please cite:
```
@misc{giammatteo2023extension,
      title={An extension of the approximate component mode synthesis method to the heterogeneous Helmholtz equation}, 
      author={Elena Giammatteo and Alexander Heinlein and Matthias Schlottbom},
      year={2023},
      eprint={2303.06671},
      archivePrefix={arXiv},
      primaryClass={math.NA}
}
```
