# Continuous Multi-Node Truss Element

This is a repository for the source codes associated with the implementation of the *continuous multi-node truss element formulation* in [OpenSees](https://github.com/OpenSees/OpenSees). 

The continuous multi-node truss element formulation enables simulation of a continuous cable/rope/wire with angular deviations while maintaining a uniform axial strain/force over its entire length. To this end, this element mormulation allows discretizing the cable/rope/wire into multiple truss sub-elements with arbitrary directions, though all these sub-elements are represented via a single multi-node element formulation and have idential axial strain/force. Both two- and three-dimensional versions of the formulation have been implemented and both use *co-rotational* geometric transformations to capture through-analysis changes to the length/direction of the truss sub-elements.

## Tcl Command

Command:
    
    element multiNodeTruss $eleTag? $nodeTags? $matTag? $A?

Arguments:

*  $eleTag    unique element object tag (integer)
*  $nodeTags  node tags from start to end (integer)
*  $matTag    uniaxial material tag (integer)
*  $A         cross-sectional area (double)

Example:

The following example constructs a multi-node truss element of tag **1** with three sub-elements connecting nodes **2**, **3**, **4**, and **5** with a constant area of **2.5** and a uniaxial material tag of **6**:

    element multiNodeTruss 1 2 3 4 5 2.5 6

## Authors

Codes written and maintained by [Mohammad Salehi](https://resilient-structures.com/) (Rice University) and [Petros Sideris](https://sites.google.com/view/petros-sideris-sem-group/) (Texas A&M University).

## References

*  Salehi, M., Sideris, P., Liel, A. (2020). "Effect of Major Design Parameters on the Seismic Performance of Bridges with Hybrid Sliding–Rocking Columns." *ASCE Journal of Bridge Engineering*, 25(10): 04020072. (https://doi.org/10.1061/(ASCE)BE.1943-5592.0001616)
*  Salehi, M. (2020). "Nonlinear Modeling, Dynamic Analysis, and Experimental Testing of Hybrid Sliding-Rocking Bridges." Ph.D. Dissertation, Zachry Department of Civil and Environmental Engineering, Texas A&M University, College Station, TX, USA. (https://hdl.handle.net/1969.1/191848)
