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

    D_{\alpha\alpha'}(q) = \braket{A_\alpha(q) A_{\alpha'}(-q)} = \frac{2\pi c^2}{\omega_{\mathbf q}}
    \Braket{\int dz \bm{\epsilon}^\ast_\alpha(\mathbf{q}, z)
        \sum_{\beta}\left(\bm{\epsilon}_{\beta}(q, z) a_{\beta,q} + \bm{\epsilon}^\ast_{\beta}(-q, z) \bar{a}_{\beta, -q}\right)
    \int dz' \bm{\epsilon}^\ast_{\alpha}(-\mathbf{q}, z')
        \sum_{\beta'}\left(\bm{\epsilon}_{\beta'}(-q, z') a_{\beta',-q} + \bm{\epsilon}^\ast_{\beta'}(q, z') \bar{a}_{\beta', q}\right)}\\
    =
    \frac{2\pi c^2}{\omega_{\mathbf q}}\int dz \epsilon^{i\ast}_\alpha(\mathbf{q}, z)\int dz' \epsilon^{j\ast}_{\alpha}(-\mathbf{q}, z')
    \sum_{\beta\beta'}
    \left(\epsilon^i_\beta(q, z) \epsilon^{\ast j}_{\beta'}(q, z')\braket{a_{\beta, q}\bar{a}_{\beta', q}})
    + \epsilon^{\ast i}_\beta(-q, z) \epsilon^{j}_{\beta'}(-q, z')\braket{\bar{a}_{\beta, -q}a_{\beta', q}})
    \right)

In Jon's notes he defines the polarizations

.. math::

    \bm{\epsilon}_1 &= \sqrt{\frac{2}{L}} \sin\left(\frac{n \pi z}{L}\right)\hat{\bm{z}} \times \hat{\bm{q}}\\
    \bm{\epsilon}_2 &= \sqrt{\frac{2}{L}} \left(\cos\left(\frac{n \pi z}{L}\right)\hat{\bm{z}} 
    - i \frac{\omega_0}{\omega_\mathbf{q}} \sin\left(\frac{n \pi z}{L}\right) \hat{\bm{q}}\right)

However, the second polarization is not properly normalized. In order to do so we must divide it by :math:`\sqrt{1 + \tfrac{\omega_0^2}{\omega_q^2}}`.
For convenience we also multiple the first polarization by i.
With thse definitions we have

.. math::

    \bm{\epsilon}_1 &= i\sqrt{\frac{2}{L}} \sin\left(\frac{n \pi z}{L}\right)\hat{\bm{z}} \times \hat{\bm{q}}\\
    \bm{\epsilon}_2 &= \sqrt{\frac{2}{L(\omega_q^2 + \omega_0^2)}} \left(\omega_\mathbf{q}\cos\left(\frac{n \pi z}{L}\right)\hat{\bm{z}} 
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

    \mathbf{A}_q(L/2) = iA_1(q)\sqrt{\frac{2}{L}} \hat{\bm{z}}\times \hat{\bm{q}} + i A_2(q)\sqrt{\frac{2}{L(\omega_q^2 + \omega_0^2)}} \omega_0 \hat{\bm{q}}
    = i \sqrt{\frac{2}{L}}\left[
        A_1(q) \hat{\bm{z}}\times \hat{\bm{q}} + A_2(q)\frac{\omega_0}{\sqrt{\omega_q^2 + \omega_0^2}} \hat{\bm{q}}\right]

This the in plane electrons couple to the two polarizations with different strengths.
This prevents us from performing a unitary transformation into a different basis in the plane.
In other words, if we wish to represent :math:`A` as components along different axes, the photonic sector will become non diagonal.

.. autodoxygenfile:: cavity.h
