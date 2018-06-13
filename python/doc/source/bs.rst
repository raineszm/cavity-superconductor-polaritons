BS
===

Mode Operators
--------------

The BS Matsubara action can be written

.. math::

    S = \frac{1}{\beta}\sum_q \left[
        \underbrace{\frac{1}{g_d} - \frac{1}{g_s}}_{m_0} - (i\Omega)^2_\text{BS}I_0((i\Omega_\text{BS})^2) \right]

where (c.f. :cpp:func:`BS::inf_gf_int()`) defines the integral :math:`I(\Omega^2)|_{v_s=0} = -\Omega^2
I_0(\Omega^2)`.

Here we have already assumed that the BS mode dispersion is irrelevant for
our purposes and noted that the supercurrent doesn't much affect the BS
mode frequency.

The condition for the BS mode frequency is then

.. math::

    m_0 - \Omega_\text{BS}^2 I_0(\Omega_\text{BS}^2) = 0

This allows us to write

.. math::

   m_0 = \Omega_\text{BS}^2 I_0(\Omega_\text{BS}^2)

We then taylor expand about the BS frequency

.. math::

    S = \frac{1}{\beta}\sum_q \left[
         \Omega_\text{BS}^2 I_0(\Omega_\text{BS}^2) - (i\Omega)^2_\text{BS}\left(I_0(\Omega_\text{BS}^2) + \left\{(i\Omega)^2
        - \Omega_\text{BS}^2\right\} I'_0(\Omega_\text{BS}^2) \right)\right]

If we only concern ourselves with low lying excitations we can throw out the derivative term and obtain
the harmonic action

.. math::

    S = \frac{1}{\beta}\sum_q \left[
         \Omega_\text{BS}^2 I_0(\Omega_\text{BS}^2) - (i\Omega)^2_\text{BS}I_0(\Omega_\text{BS}^2)\right]

Defining :math:`I_0(\Omega_\text{BS}^2) = \frac{1}{2}M` we have a Harmonic oscilator and we can immediately write

.. math::

    d_q = \frac{b_q + b^\dagger_{-q}}{\sqrt{2 M \Omega_\text{bs}}}
    = \frac{b_q + b^\dagger_{-q}}{\sqrt{4 I_0(\Omega_\text{BS}) \Omega_\text{BS}}}
    = \frac{1}{2}\sqrt{\frac{\Omega_\text{BS}}{m_0}}\left(b_q + b^\dagger_{-q}\right)

.. autodoxygenfile:: bs.h
