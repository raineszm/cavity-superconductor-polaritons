Coupling
=========

Coupling in the Mode Basis
============================

As detailed in :doc:`cavity`, in terms of the mode operators

.. math::

    \mathbf{A}(\mathbf{q}, z) = \sqrt{\frac{2\pi c^2}{\omega_\mathbf{q}}}\left(\bm{\epsilon}_\alpha(q, z) a_{\alpha,q} + \bm{\epsilon}_\alpha(-q, z)^* a^\dagger_{-q, \alpha}\right)

and the mode action is

.. math::

    S_a = \sum_q \bar{a}_{\alpha, q}(-i \omega_m + \omega_\mathbf{q}) a_{\alpha, q}

Consider then the coupling of the BS mode to the photons in terms of mode operators (c.f. :cpp:func:`BS::hamiltonian`, :cpp:func:`Coupling::ImDA`)

This coupling can be written

 .. math::

    \sum_q i g(i \Omega, \mathbf{q}) \left(d_{-q} \mathbf{v}_s \cdot \mathbf{A}_q(L/2) - \mathbf{A}_{-q}(L/2)\cdot \mathbf{v}_s d_q\right)

where :math:`g(\Omega, \mathbf{q})` is the coupling given by :cpp:func:`Coupling::ImDA`.
We also make use of the BS mode operators (c.f. :doc:`bs`)

.. math::

    d_q = \frac{b_q + b^\dagger_{-q}}{\sqrt{2 M \Omega_\text{bs}}}
    = \frac{b_q + b^\dagger_{-q}}{\sqrt{4 I_0(\Omega_\text{BS}) \Omega_\text{BS}}}
    = \frac{1}{2}\sqrt{\frac{\Omega_\text{BS}}{m_0}}\left(b_q + b^\dagger_{-q}\right)


Substituting into the coupling term

 .. math::

   2i\sqrt{\frac{\Omega_\text{BS} \pi c^2}{2 m_0}} \sum_q i g(i \Omega, \mathbf{q})\frac{1}{\sqrt{\omega_\mathbf{q}}}
   \sum_\alpha
   \left[
     \left(b_{-q} + b^\dagger_{q}\right)\mathbf{v}_s \cdot \left(\bm{\epsilon}_\alpha(q, L/2) a_{\alpha,q} + \bm{\epsilon}_\alpha(-q, L/2)^* a^\dagger_{-q, \alpha}\right)
   \right]

where we have used the symmetry of the integrand with respect to :math:`q\to-q`


Throwing out counter-rotating terms


 .. math::

   2i\sqrt{\frac{\Omega_\text{BS} \pi c^2}{2 m_0}} \sum_q i g(i \Omega, \mathbf{q})\frac{1}{\sqrt{\omega_\mathbf{q}}}
   \sum_\alpha
   \left[
     b^\dagger_{q}\mathbf{v}_s \cdot \bm{\epsilon}_\alpha(q, L/2) a_{\alpha,q}
     + a^\dagger_{-q, \alpha}\mathbf{v}_s\bm{\epsilon}_\alpha(-q, L/2)^*b_{-q}
   \right]

Using the properties of :math:`\bm{\epsilon}` this can be written


.. math::
   \sqrt{\frac{2\pi c^2\Omega_\text{BS} }{m_0}} \sum_q  g(i \Omega, \mathbf{q})\frac{1}{\sqrt{\omega_\mathbf{q}}}
   \sum_\alpha
   \left[
     b^\dagger_{q}\mathbf{v}_s \cdot i\bm{\epsilon}_\alpha(q, L/2) a_{\alpha,q}
     + h.c.
   \right]

We may then extract the effective coupling

.. math::

    g_\text{eff}(q, \alpha) = \sqrt{\frac{2\pi c^2\Omega_\text{BS} }{m_0\omega_{\mathbf{q}}}}
    \mathbf{v}_s\cdot i\bm{\epsilon}_{\alpha,\mathbf{q}}\left(\frac{L}{2}\right)g(i\Omega, \mathbf{q})

or explicitly

.. math::


   g_\text{eff}(i\Omega_m, q) \approx -2i \Omega_m e \mathbf{v}_s\cdot i\bm{\epsilon}_{\alpha,\mathbf{q}}
   \nu \Delta \sqrt{\frac{2\pi \Omega_\text{BS} }{m_0\omega_{\mathbf{q}}}}
    \int_\Delta^\infty
   \frac{d\lambda}{\sqrt{\lambda^2 - \Delta^2}}
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
   \tilde{\Pi}_{\alpha\beta}(q) = \frac{1}{2}\frac{\pi c^2}{\omega_q}\bm{\epsilon}_\alpha(q, L/2)^*
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


.. autodoxygenfile:: coupling.h
