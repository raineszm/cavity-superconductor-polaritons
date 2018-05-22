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

    \sqrt{\frac{2}{L}}\mathbf{v} \cdot \begin{pmatrix}A_x(q)\\A_y(q)\end{pmatrix}

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

.. autodoxygenfile:: cavity.h
