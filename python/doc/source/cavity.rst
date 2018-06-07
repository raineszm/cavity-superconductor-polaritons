Cavity
=======


Vector Potential and Mode Operators
--------------------------------------

In this code we use Gaussian (unrationalized) units.

In the cavity our field can be expressed in terms of mode operators as

.. math::

    \mathbf{A}(\mathbf{q}, z) = C\left(\mathbf{\epsilon}_\alpha(q, z) a_{\alpha,q} + \mathbf{\epsilon}_\alpha(-q, z)^* a^\dagger_{-q, \alpha}\right)

It will however be useful to define the fields

.. math::

    \mathbf{A}_\alpha(\mathbf{q}) = C_\mathbf{q}\int_0^L dz \bm{\epsilon}^\ast_\alpha(\mathbf{q}, z)
    \sum_{\alpha'}\left(\bm{\epsilon}_{\alpha'}(q, z) a_{\alpha',q} + \bm{\epsilon}^\ast_{\alpha'}(-q, z) a^\dagger_{-q, \alpha'}\right)

The prefactor :math:`C_\mathbf{q}` normalizes the energy so that

.. math::

    u = \frac{1}{8\pi} \int d\mathbf{x} \left(E^2 + B^2\right) = \sum_\mathbf{q} \omega_\mathbf{q} a_{q,\alpha}^\dagger a_{q, \alpha}

where :math:`\omega^2_\mathbf{q}= \omega_0^2 + c^2 q^2`.

In Gaussian units and Coulomb gauge (with :math:`\phi=0`)

.. math::

    \mathbf{E} = -\frac{1}{c}\frac{\partial}{\partial t} \mathbf{A}\\
    \mathbf{B} = \nabla \times \mathbf{A}

We can obtain :math:`\mathbf{E}` from the Heisenberg equation of motion

.. math::

    d_t \mathbf{A} = -c \mathbf{E} = i [H, \mathbf{A}]

Using the relation :math:`[a_{q, \alpha}, a^\dagger_{q', \alpha'}] = (2\pi)^2 \delta(q-q') \delta_{\alpha\alpha'}`

.. math::

    \mathbf{E}(\mathbf q, z) = - \frac{i}{c} C \omega_q \left( \bm{\epsilon}_\alpha(q, z)  a_{\mathbf q, \alpha} - \bm{\epsilon}_\alpha(-q, z)^* a^\dagger_{-\mathbf{q},\alpha}\right)

So

.. math::

    \int d\mathbf{x} \left| E \right|^2 =  2\sum_q  \frac{C_\mathbf{q}^2}{c^2} \omega_q^2 a^\dagger_{q, \alpha} a_{q, \alpha}

up to ordering ambiguity

For a free space solution, the magnetic and electric energies will be equal so

.. math::

    u = \frac{1}{2\pi c^2} \sum_q C_\mathbf{q}^2 \omega_q^2 a^\dagger_{q, \alpha} a_{q, \alpha} = \sum _q \omega_q a^\dagger_{q, \alpha} a_{q, \alpha}
    \implies C_\mathbf{q}^2 = \frac{2\pi c^2}{\omega_\mathbf{q}}


So the relation between mode operators and the A field is

.. math::

    \mathbf{A}(\mathbf{q}, z) = \sqrt{\frac{2\pi c^2}{\omega_\mathbf{q}}}\left(\bm{\epsilon}_\alpha(q, z) a_{\alpha,q} + \bm{\epsilon}_\alpha(-q, z)^* a^\dagger_{-q, \alpha}\right)


Thermal Photon Action
----------------------

Beginning with the mode action we have

.. math::

    S_a = \sum_q \bar{a}_{\alpha, q}(-i \omega_m + \omega_\mathbf{q}) a_{\alpha, q}

On the other hand we have an equivalent description in terms of vector potential

.. math::

    S_A = \frac{1}{2}\sum_q A_\alpha(-q) D^{-1}_{\alpha\alpha'}(q) A_{\alpha'}(q)

where :math:`D^{-1}` in the inverse photon propagator. We can relate the two actions by

.. math::
    :nowrap:

    \begin{multline}
    D_{\alpha\alpha'}(q) = \braket{A_\alpha(q) A_{\alpha'}(-q)}\\
     = \frac{2\pi c^2}{\omega_{\mathbf q}}
    \left\langle\int dz \bm{\epsilon}^\ast_\alpha(\mathbf{q}, z)
        \sum_{\beta}\left(\bm{\epsilon}_{\beta}(q, z) a_{\beta,q} + \bm{\epsilon}^\ast_{\beta}(-q, z) \bar{a}_{\beta, -q}\right)\right.\\
        \times\left.
    \int dz' \bm{\epsilon}^\ast_{\alpha}(-\mathbf{q}, z')
        \sum_{\beta'}\left(\bm{\epsilon}_{\beta'}(-q, z') a_{\beta',-q} + \bm{\epsilon}^\ast_{\beta'}(q, z') \bar{a}_{\beta', q}\right)\right\rangle\\
    =
    \frac{2\pi c^2}{\omega_{\mathbf q}}\int dz \epsilon^{i\ast}_\alpha(\mathbf{q}, z)\int dz' \epsilon^{j\ast}_{\alpha}(-\mathbf{q}, z')
    \sum_{\beta\beta'}
    \left(\epsilon^i_\beta(q, z) \epsilon^{\ast j}_{\beta'}(q, z')\braket{a_{\beta, q}\bar{a}_{\beta', q}})
    + \epsilon^{\ast i}_\beta(-q, z) \epsilon^{j}_{\beta'}(-q, z')\braket{\bar{a}_{\beta, -q}a_{\beta', q}})
    \right)
    \end{multline}

The n-th harmonic TE and TM polarizations are respectively

.. math::

    \bm{\epsilon}_1 &= i\sqrt{\frac{2}{L}} \sin\left(\frac{n \pi z}{L}\right)\hat{\bm{z}} \times \hat{\bm{q}}\\
    \bm{\epsilon}_2 &= \sqrt{\frac{2}{L}} \frac{1}{\omega_\mathbf{q}}\left(c q\cos\left(\frac{n \pi z}{L}\right)\hat{\bm{z}}
    -i \omega_0 \sin\left(\frac{n \pi z}{L}\right) \hat{\bm{q}}\right)

Note that both of these have the property :math:`\bm{\epsilon(-q, z)} = \bm{\epsilon}^\ast(q, z)`.
And from the mode action

.. math::

    \braket{a_{\beta, q}\bar{a}_{\beta', q}}) = \frac{1}{i \omega_m - \omega_\mathbf{q}} \delta_{\beta\beta'}

This allows us to write

.. math::


    D_{\alpha\alpha'}(q)
    = \frac{2\pi c^2}{\omega_{\mathbf q}}\int dz \epsilon^{i\ast}_\alpha(\mathbf{q}, z)\int dz' \epsilon^{j}_{\alpha}(\mathbf{q}, z')
    \sum_{\beta}
    \left(\epsilon^i_\beta(q, z) \epsilon^{\ast j}_{\beta}(q, z') \frac{1}{i\omega_m - \omega_\mathbf{q}}
    + \epsilon^{i}_\beta(q, z) \epsilon^{j\ast}_{\beta}(q, z')\frac{1}{-i\omega_m - \omega_\mathbf{q}}
    \right)

Using the orthonormality of the polarizations this becomes


.. math::

    D_{\alpha\alpha'}(q) =  \frac{2\pi c^2}{\omega_{\mathbf q}} \frac{2 \omega_\mathbf{q}}{(i\omega_m)^2 - \omega_\mathbf{q}^2} \delta_{\alpha\alpha'}

We may then immediately invert this to obtain the action

.. math::

    S_A = \frac{1}{8 \pi c^2}\sum_q A_\alpha(-q) \left[ (i \omega_m)^2 - \omega_\mathbf{q}^2\right]A_{\alpha'}(q)

The question then remains how these :math:`A_\alpha` fields couple to fermions.

We recall

.. math::

    \mathbf{A}_\alpha(\mathbf{q}) = \sqrt{\frac{2\pi c^2}{\omega_\mathbf{q}}}
    \left(a_{\alpha,q} + a^\dagger_{-q, \alpha}\right)

and

.. math::

    \mathbf{A}(\mathbf{q}, z) = \sqrt{\frac{2\pi c^2}{\omega_\mathbf{q}}}\sum_\alpha\bm{\epsilon}_\alpha(q, z) \left(a_{\alpha,q} + a^\dagger_{-q, \alpha}\right)

where we have used the transformation properties of :math:`\epsilon` under :math:`q \to -q`.
We immediately see the vector potential can be expressed as :math:`\mathbf{A} = \sum_\alpha \bm{\epsilon_\alpha} A_\alpha`.
In the plane we then have

.. math::

    \mathbf{A}_q(L/2)
    = i \sqrt{\frac{2}{L}}\left[
        A_1(q) \hat{\bm{z}}\times \hat{\bm{q}} - A_2(q)\frac{\omega_0}{\omega_q} \hat{\bm{q}}\right]

Thus the in plane electrons couple to the two polarizations with different strengths.
This prevents us from performing a unitary transformation into a different basis in the plane.
In other words, if we wish to represent :math:`A` as components along different axes, the photonic sector will become non diagonal.

Suppose we wish to write the theory such that the paramagnetic coupling is


.. math::

    \frac{e}{c}\sqrt{\frac{2}{L}}\mathbf{v} \cdot \begin{pmatrix}A_x(q)\\A_y(q)\end{pmatrix}

We can do this by defining the transformation

.. math::

    \begin{pmatrix}
    A_x
    A_y
    \end{pmatrix} = \underbrace{-i
    \begin{pmatrix}
    \sin \theta_q& \frac{\omega_0}{\omega_\mathbf{q}} \cos \theta_q \\
    -\cos \theta_q&  \frac{\omega_0}{\omega_\mathbf{q}} \sin\theta_q
    \end{pmatrix}}_{U(q)}
    \begin{pmatrix}
    A_1(q)\\
    A_2(q)
    \end{pmatrix}

This is a non-unitary transformation but since it is linear the contribution to the Jacobian cancels out in any expectation value.
However, the cavity action in this basis becomes

.. math::

    S_A = \frac{1}{8 \pi c^2}\sum_q \mathbf{A}(-q) \left[ (i \omega_m)^2 - \omega_\mathbf{q}^2\right](U^{-1}(-\mathbf q))^T U^{-1}(\mathbf q)\mathbf{A}(q)\\
    = \frac{1}{8 \pi c^2}\sum_q \mathbf{A}(-q) \left[ (i \omega_m)^2 - \omega_\mathbf{q}^2\right]
    \begin{pmatrix}
    \left(\frac{\omega_q}{\omega_0}\right)^2 \cos^2 \theta_q + \sin^2 \theta_q& \left(\frac{\omega_q^2}{\omega_0^2} -1\right)\sin\theta_q \cos\theta_q\\
    \left(\frac{\omega_q^2}{\omega_0^2} -1\right)\sin\theta_q \cos\theta_q&\left(\frac{\omega_q}{\omega_0}\right)^2 \sin^2 \theta_q + \cos^2 \theta_q
    \end{pmatrix}
    \mathbf{A}(q)

Finally if we wish to transform to the basis along and perpendicular to the supercurrent

.. math::

    \begin{pmatrix}
    A^x\\
    A^y\\
    \end{pmatrix}
    = \underbrace{\begin{pmatrix}
    \cos \theta_s& -\sin\theta_s\\
    \sin\theta_s& \cos\theta_s
    \end{pmatrix}}_R
    \begin{pmatrix}
    A_\parallel\\
    A_\perp
    \end{pmatrix}

Upon this transformation the photon action becomes

.. math::

    S_A = \frac{1}{16 \pi c^2}\sum_q \mathbf{A}(-q) \left[ (i \omega_m)^2 - \omega_\mathbf{q}^2\right]
    \left[
    \left(1 + \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
    - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin 2(\theta_q - \theta_s)\sigma_1 - \cos 2(\theta_q - \theta_s)\sigma_3\right)
    \right]
    \mathbf{A}(q)

On the other hand.
We generally want all terms in the inverse propagator to have the same units.
Looking at the coupling term (c.f. :cpp:func:`Coupling::ImDA`)

.. math::

    g_q \approx -2 \sqrt{\frac{2}{L}}\frac{e}{c} \nu v_s \Omega\Delta \int_\Delta^\infty
   \frac{d\lambda}{\sqrt{\lambda^2 - \Delta^2}}
   \int_0^{2\pi}\frac{d\theta}{2\pi}
   \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(\Omega + i0^+)^2 -
   (2\lambda)^2}f_d(\theta)

By inspection the units of this term are

.. math::
    [g_q] = [e \nu \frac{v_s}{c}/\sqrt{L}]

By inspecting the BS action (c.f. :cpp:func:`BS::action`) one can see that the Bardasis-Schrieffer inverse GF
has the same units as :math:`\nu`.
As such it makes sense to absorb the factor :math:`\sqrt{\tfrac{2}{L}}e` into the photon fields.
This makes the paramagnetic coupling

.. math::

    \frac{\mathbf{v}}{c} \cdot \mathbf{A}

and the coupling between the excitations

.. math::

    g_q \approx -2 \frac{v_s}{c} \nu \Omega\Delta \int_\Delta^\infty
   \frac{d\lambda}{\sqrt{\lambda^2 - \Delta^2}}
   \int_0^{2\pi}\frac{d\theta}{2\pi}
   \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(\Omega + i0^+)^2 -
   (2\lambda)^2}f_d(\theta)

and the photon action


.. math::

    S_A = \frac{L}{32\pi e^2 c^2}\sum_q \mathbf{A}(-q) \left[ (i \omega_m)^2 - \omega_\mathbf{q}^2\right]
    \left[
    \left(1 + \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
    - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin 2(\theta_q - \theta_s)\sigma_1 - \cos 2(\theta_q - \theta_s)\sigma_3\right)
    \right]
    \mathbf{A}(q)

or using :math:`\alpha=\frac{e^2}{c}`

.. math::

    S_A = \frac{\alpha^2 L}{32\pi (\alpha c)^3}\sum_q \mathbf{A}(-q) \left[ (i \omega_m)^2 - \omega_\mathbf{q}^2\right]
    \left[
    \left(1 + \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
    - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin 2(\theta_q - \theta_s)\sigma_1 - \cos 2(\theta_q - \theta_s)\sigma_3\right)
    \right]
    \mathbf{A}(q)

Returning to the mode Basis
============================

As detailed above, in terms of the mode operators

.. math::

    \mathbf{A}(\mathbf{q}, z) = \sqrt{\frac{2\pi c^2}{\omega_\mathbf{q}}}\left(\bm{\epsilon}_\alpha(q, z) a_{\alpha,q} + \bm{\epsilon}_\alpha(-q, z)^* a^\dagger_{-q, \alpha}\right)

and the mode action is

.. math::

    S_a = \sum_q \bar{a}_{\alpha, q}(-i \omega_m + \omega_\mathbf{q}) a_{\alpha, q}

Consider then the coupling of the BS mode to the photons in terms of mode operators (c.f. :cpp:func:`BS::hamiltonian`, :cpp:func:`Coupling::ImDA`)

 This coupling can be written

 .. math::

    \sum_q i g(q) \left(d_{-q} \mathbf{v}_s \cdot \mathbf{A}_q(L/2) - \mathbf{A}_{-q}(L/2)\cdot \mathbf{v}_s d_q\right)

We now substitute in mode operators

 .. math::

    i\sqrt{\frac{2\pi c^2}{2M\Omega_\text{BS}}}
    \sum_q \sqrt{\frac{1}{\omega_q}}g(q) \\
    \times
    \left[\left(b_{-q} + b^\dagger_{q}\right)\mathbf{v}_s \cdot \sum_\alpha \bm{\epsilon}_{q,\alpha}(L/2)\left(a_{\alpha,q} + a^\dagger_{\alpha, -q}\right)\right.\\
    \left.-\mathbf{v}_s \cdot \sum_\alpha \bm{\epsilon}_{q,\alpha}(L/2)\left(a_{\alpha,-q} + a^\dagger_{\alpha,q}\right)\left(b_q + b^\dagger_{-q}\right)\right]

Throwing out counter-rotating terms we can write this


 .. math::

    i\sqrt{\frac{\pi c^2}{M\Omega_\text{BS}}}
    \sum_q \sqrt{\frac{1}{\omega_q}}g(q)
    \mathbf{v}_s \cdot \sum_\alpha \bm{\epsilon}_{q,\alpha}(L/2)
    \left[a^\dagger_{\alpha, -q} b_{-q} + b^\dagger_{q}a_{\alpha,q}
    -b^\dagger_{-q}a_{\alpha,-q} - a^\dagger_{\alpha,q}b_q\right]

And using the symmetry properties of :math:`\epsilon`

 .. math::

    S_g = 2i\sqrt{\frac{\pi c^2}{M\Omega_\text{BS}}}
    \sum_q \frac{1}{\sqrt{\omega_q}}g(q)
    \mathbf{v}_s \cdot \sum_\alpha \bm{\epsilon}_{q,\alpha}(L/2)
    \left[b^\dagger_{q}a_{\alpha,q} -a^\dagger_{\alpha,q}b_q\right]

We may then extract the effective coupling

.. math::

    g_\text{eff}(q, \alpha) = 2 i \sqrt{\frac{\pi c^2}{M\Omega_\text{BS}\omega_q}}
    \mathbf{v}_s\cdot \bm{\epsilon}_{\alpha,\mathbf{q}}\left(\frac{L}{2}\right)g(q)

or explicitly

.. math::

    g_\text{eff}(q, \alpha) = 4 i e\nu\Omega\Delta  \sqrt{\frac{\pi}{M\Omega_\text{BS}\omega_q}}
    \mathbf{v}_s\cdot \bm{\epsilon}_{\alpha,\mathbf{q}}\left(\frac{L}{2}\right)
    \int_\Delta^\infty
    \frac{d\lambda}{\sqrt{\lambda^2 - \Delta^2}}\\
    \times
    \int_0^{2\pi}\frac{d\theta}{2\pi}
    \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(\Omega + i0^+)^2 -
    (2\lambda)^2}f_d(\theta)


Photon Self-Energy
==================

Consider the addition of the photonic self-energy to the photon action.

The thermal photon action is (c.f. :cpp:func:`Coupling::photon_se`)

.. math::

    \mathcal{L}_\text{SE} = -\frac{1}{2\beta}\sum_q\mathbf{A}_{-q}(L/2)\hat{\Pi}(q)\mathcal{A}_q(L/2)

Going to the mode basis

.. math::

    \frac{1}{2}\sum_{\alpha,\beta}
    \frac{\pi c^2}{\omega_\mathbf{q}}
    \left(\bm{\epsilon}_\alpha(-q, L/2) a_{\alpha,-q} + \bm{\epsilon}_\alpha(q, L/2)^* a^\dagger_{q, \alpha}\right)
    \hat{\Pi}(q)
    \left(\bm{\epsilon}_\beta(q, L/2) a_{\beta,q} + \bm{\epsilon}_\beta(-q, L/2)^* a^\dagger_{-q, \beta}\right)

Ignoring the counter rotating terms

.. math::

    \frac{1}{2}\sum_{\alpha,\beta}
    \frac{\pi c^2}{\omega_\mathbf{q}}\left(
    a^\dagger_{q, \alpha}\bm{\epsilon}_\alpha(q, L/2)^*
    \hat{\Pi}(q)
    \bm{\epsilon}_\beta(q, L/2) a_{q,\beta}
    +
    a_{\alpha,-q}\bm{\epsilon}_\alpha(-q, L/2)
    \hat{\Pi}(q)
    \bm{\epsilon}_\beta(-q, L/2)^* a^\dagger_{-q, \beta }
    \right)\\ =
    \frac{1}{2}\sum_{\alpha,\beta}
    \frac{\pi c^2}{\omega_\mathbf{q}}
    a^\dagger_{q, \alpha}\bm{\epsilon}_\alpha(q, L/2)^*
    \left(
    \hat{\Pi}(q)
    +
    \hat{\Pi}^T(-q)\right)
    \bm{\epsilon}_\beta(q, L/2)
    a_{q,\beta}

We thus define

.. math::
   \tilde{\Pi}_{\alpha\beta}(q) = \frac{1}{2}\bm{\epsilon}_\alpha(q, L/2)^*
    \left(
    \hat{\Pi}(q)
    +
    \hat{\Pi}^T(-q)\right)
    \bm{\epsilon}_\beta(q, L/2)

The thermal photon action is then

.. math::

    S = \frac{1}{\beta}\sum_q a^\dagger_{q, \alpha}\left(-i\omega_m + \omega_\mathbf{q} + \tilde{\Pi}_{\alpha\beta}(i \omega_m, \mathbf{q})\right)a_{q,\beta}

Renormalization
---------------

In order to normalize we must first find the new mass.
At :math:`q=0`

.. math::

   S = -i\omega_m + \omega_0 + \tilde{\Pi}_{\alpha\beta}(i \omega_m, 0)

The renormalized mass :math:`\omega_r` is the frequency at which this action vanishes.
This allows us to expand

.. math::

    \tilde{\Pi} \approx (\hat{Z}-1)\left(\omega_r - i \omega_m  - \omega_0\right) + \hat{\tilde{\Pi}}(\omega_r, \mathbf{q}) + \cdots

where

.. math::

    1 - \hat{Z} = \left.\frac{\partial\Pi(i\omega, 0)}{\partial(i\omega)}\right|_{i\omega=\omega_r-\omega_0}


Assuming :math:`\hat{Z}` is positive definite it admits a Cholesky decomposition :math:`\hat{Z} = \hat{L} \hat{L}^\dagger`.
We then absorb the matrix :math:`\hat{L}` in the definition of our field operators

.. math::

    \tilde{a} = \hat{L}^\dagger a

This makes the photonic Lagrangian

.. math::

    \mathcal{L} = -i\omega_m + \tilde{\omega_\mathbf{q}} + \hat{L}^{-1}\hat{\tilde{\Pi}}(\omega_r, \mathbf{q})(\hat{L}^{\dagger})^{-1}

allowing us to define the effective Hamiltonian

.. math::

    \hat{H}_\text{phot} = \tilde{\omega}_{\mathbf{q}} + \hat{L}^{-1}\hat{\tilde{\Pi}}(\omega_r, \mathbf{q}){(\hat{L}^\dagger)}^{-1}

Similarly, the coupling to the Bardasis-Schrieffer mode becomes

.. math::

    \tilde{g}_\text{eff}(q, \alpha) = 4 i e\nu\Omega\Delta  \sqrt{\frac{\pi}{M\Omega_\text{BS}\omega_q}}
    \mathbf{v}_s\cdot \bm{\epsilon}_{\alpha',\mathbf{q}} \left[{(\hat{L}^\dagger)}^{-1}\right]_{\alpha'\alpha} \left(\frac{L}{2}\right)
    \int_\Delta^\infty
    \frac{d\lambda}{\sqrt{\lambda^2 - \Delta^2}}\\
    \times
    \int_0^{2\pi}\frac{d\theta}{2\pi}
    \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(\Omega + i0^+)^2 -
    (2\lambda)^2}f_d(\theta)


.. autodoxygenfile:: cavity.h

