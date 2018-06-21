BS
===

Mode Operators
--------------

The BS Matsubara action can be written

.. math::

    S = \frac{1}{\beta}\sum_q \left[
        \underbrace{\frac{1}{g_d} - \frac{1}{g_s}}_{m_0} + I((i\Omega_\text{BS})^2) \right]

where (c.f. :cpp:func:`BS::inv_gf_int()`) defines the integral :math:`I(\Omega^2)`.

Here we have already assumed that the BS mode dispersion is irrelevant for
our purposes and noted that the supercurrent doesn't much affect the BS
mode frequency.

The condition for the BS mode frequency is then

.. math::

    m_0 +I(\Omega_\text{BS}^2) = 0

This allows us to write

.. math::

   m_0 = -I(\Omega_\text{BS}^2)

We then taylor expand about the BS frequency

.. math::

    S = \frac{1}{\beta}\sum_q d_{-q}\left[
        \underbrace{m_0 + I(\Omega_\text{BS}^2)}_0 + \left[(i\Omega_m)^2 - \Omega_\text{BS}^2\right]I'(\Omega_\text{BS}^2) + \cdots\right]d_q

If we only concern ourselves with low lying excitations we obtain
the harmonic action

.. math::

    S = \frac{1}{\beta}\sum_q 
    d_{-q} \left[(i\Omega_m)^2 - \Omega_\text{BS}^2\right]I'(\Omega_\text{BS}^2)d_q

Defining :math:`I'(\Omega_\text{BS}^2) = -\frac{1}{2}M` we have a Harmonic oscilator and we can immediately write

.. math::

    d_q = \frac{b_q + b^\dagger_{-q}}{\sqrt{2 M \Omega_\text{bs}}}
    = \frac{b_q + b^\dagger_{-q}}{\sqrt{4 I_0(\Omega_\text{BS}) \Omega_\text{BS}}}
    = \frac{1}{2}\sqrt{\frac{\Omega_\text{BS}}{m_0}}\left(b_q + b^\dagger_{-q}\right)

.. autodoxygenfile:: bs.h
