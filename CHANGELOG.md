# Changelog

## 0.2.0
- Spin correlations are now measured with respect to all lattice basis site within a single unit cell. 
This allows for smoother calculation of structure factors (in particular if not all lattice sites are symmetry equivalent). See also issue #7. 
The previous behavior, which computed correlations only with respect to the first basis site, can be reproduced by accessing only the first column of the matrix. For example, for a MonteCarlo simulation instance `mc`, the mean correlation with respect to the first basis site could be obtained as `mean(mc.observables.correlation)[:,1]`.
- Fixed a missing check that would allow instantiating a lattice with self-interactions, which leads to incorrect calculations of interaction energy. Such check already existed for on-site interactions, but not for self-interactions that result from crossing one or more periodic boundaries. See also issue #8.
