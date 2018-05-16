Summary {#mainpage}
=========

This code is for the paper *"Enriched axial anomaly in Weyl materials"*.
Its main purpose is to calculate the non-linear susceptibility of a lattice model to simultaneous 'charge'- and 'spin'-density waves.

The lattice model consists of a 'free' part
\f[
    H_0 = \sum_{\mathbf k} c^\dagger_\mathbf{k}
\left[
    \epsilon(\mathbf k) + \mathbf{d}(\mathbf k) \cdot \sigma
\right]
c_{\mathbf k}
\f]
with \f$\epsilon(\mathbf k) = t_1 \sin k_z\f$ and \f$\mathbf{d}(\mathbf k) = \left(\sin k_x,\sin k_y,2 + \cos b_z - \sum_i\cos k_i \right)^T\f$,
which is host to a pair of Weyl fermions. The momentum-space separation of the nodes is given by \f$2\mathbf{b} = 2b_z\mathbf{\hat z}\f$ and
the energy separation by \f$2b_0 = 2t_1 \sin b_z\f$.
Note that the $\sin k_z$ term also tilts the low energy theory to a type II WSM. To alleviate this we can cancel off that tilt by adding a term proportional to \f$\sin 2 k_z\f$. For more information see Params::adjusted_ccp and eps().

We then wish to calculate the lowest order response of the current to the perturbations
\f[
    H_m = m\sum_{\mathbf k} c^\dagger_{\mathbf k + \mathbf b}e^{-i\alpha}\sigma_z c_{\mathbf k - \mathbf b} + {\rm h.c.},
\f]
and
\f[
    H_g = \sum_{\mathbf k} c^\dagger_{\mathbf k + \mathbf b}\sigma_z \mathbf{g}(\tau)\cdot\sigma c_{\mathbf k - \mathbf b} + {\rm h.c.}.
\f]
