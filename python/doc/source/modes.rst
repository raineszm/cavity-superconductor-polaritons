Modes and Mode Operators
=========================

Consider a real field defined with the Lagrangian

.. math::

    \mathcal{L} =  \frac{1}{2} m |\dot{\phi}_q|^2 -\frac{1}{2} m \omega_q^2 |\phi_q|^2

where :math:`\phi^*_q = \phi_{-q}`.
How do we go to the mode basis in the path integral picture.
We start by introducing dual field :math:`\Pi` with Lagrangian

.. math::

    \mathcal{L}_\Pi = -\frac{1}{2m}|\Pi_q^2|

We now make the shift :math:`\Pi_q \to \Pi_q - m\dot{\phi}_q`.

.. math::

    \mathcal{L} \to - \frac{1}{2} m \omega_q^2 |\phi_q|^2 + \dot{\phi}_q \Pi_{-q} - \frac{1}{2m} |\Pi_q|^2

We now make the following change of basis

.. math::

    b_q = \frac{1}{\sqrt{2m\omega_q}}\left(m\omega_q\phi_q + i \Pi_q\right)\\
    b^\dagger_q = \frac{1}{\sqrt{2m\omega_q}}\left(m\omega_q\phi_{-q} - i \Pi_{-q}\right)

where we have assumed :math:`\omega_q` to be even in q.
This transformation has determinant 1.

Inverting the transformation

.. math::

    \phi_q = \frac{b_q + b^\dagger_{-q}}{\sqrt{2m\omega_q}}\\
    \Pi_q = i\sqrt{\frac{m \omega_q}{2}} \left(b^\dagger_{-q} - b_q\right)

And inserting this into the Lagrangian


.. math::

    \mathcal{L} \to
    - \frac{1}{2} m \omega_q^2 \frac{1}{2m\omega_q} \left(b_q + b^\dagger_{-q}\right)\left(b^\dagger_q + b_{-q}\right)\\
    - \frac{1}{2m} \frac{m \omega_q}{2}\left(b^\dagger_{-q} - b_q\right) \left(b^\dagger_{q} - b_{-q}\right)\\
    + \frac{i}{2} \left(b^\dagger_{q} - b_{-q}\right) \left(\dot{b}^\dagger_{-q} + \dot{b}_q\right)

After expanding (up to total derivatives) we have

.. math::

    \mathcal{L} = i b^\dagger_{q} \partial_t b_q - \omega_q b^\dagger_q b_q

For a Harmonic-like potential we can go to the mode basis by replacing the quadratic action for the field :math:`\phi` by
that for :math:`b` and replace :math:`\phi` everywhere else by


.. math::

    \phi_q = \frac{b_q + b^\dagger_{-q}}{\sqrt{2m\omega_q}}

In Matsubara Time
=================

Consider the Euclidean action

.. math::

    \int d\tau \sum_q \bar{b}_q (\partial_\tau + \omega_q) b_q

Let us parametrize in terms of the real and imagniary parts of :math:`b`

.. math::

    b_q = X_q + i Y_q\\
    \bar{b}_q = X_{-q} - i Y_{-q}

Plugging this in gives

.. math::

    \int d\tau \sum_q (X_{-q} - i Y_{-q}) (\partial_\tau + \omega_q)
    (X_q + i Y_q) \\
    =
    \int d\tau \sum_q \left(i X_{-q} \partial_\tau Y_q - i Y_q \partial_\tau X_{-q}
    + \omega_q (|X_q|^2 + |Y_q|^2)
    \right)

where we have ignored complete derivatives due to our boundary conditions.

Or integrating by parts

.. math::

    \int d\tau \sum_q \left(- 2 i Y_q \partial_\tau X_{-q}
    + \omega_q (|X_q|^2 + |Y_q|^2)
    \right)

In order to obtain the usual harmonic potential we define

.. math::

    X_q = \sqrt{\frac{m \omega_q}{2}}\phi_q

.. math::

    \int d\tau \sum_q \left(-  i \sqrt{2 m \omega_q}Y_q \partial_\tau \phi_{-q}
    + \frac{1}{2}m \omega^2_q |\phi_q|^2 + \omega_q|Y_q|^2
    \right)

In order to make our change of variables unitary we then define

.. math::

    Y_q = \frac{1}{\sqrt{2 m \omega_q}} \Pi_q

For

.. math::

    \int d\tau \sum_q \left(-  i \Pi_q \partial_\tau \phi_{-q}
    + \frac{1}{2}m \omega^2_q |\phi_q|^2 + \frac{1}{2m}\Pi_q\Pi_{-q}
    \right)

We now shift :math:`\Pi_q \to \Pi_q + i m \partial_\tau \phi_q` which gives

.. math::

    \int d\tau \sum_q \left(
    i (-i)m |\partial_\tau \phi_q|^2
    - \frac{1}{2}m |\partial_\tau \phi_q|^2
    + \frac{1}{2}m \omega^2_q |\phi_q|^2 + \frac{1}{2m}\Pi_q \Pi_{-q}
    \right) \\
    =
    \int d\tau \sum_q \left(
    \frac{1}{2} m |\partial_t \phi_q|^2
    + \frac{1}{2}m \omega^2_q |\phi_q|^2 + \frac{1}{2m}\Pi_q \Pi_{-q}
    \right)


.. note:: We have deformed the contour for :math:`\Pi` off of the real line. Regardless, due to the lack of poles we are free to shift the contour back down to the axis.


We are then free to integrate out :math:`\Pi_q` which gives a non-essential constant and we obtain

.. math::

    \int d\tau \sum_q \left(
    \frac{1}{2} m |\partial_\tau \phi_q|^2
    + \frac{1}{2}m \omega^2_q |\phi_q|^2
    \right)

Thus we see that in terms of the mode operators the Matsubara field is

.. math::

    \phi_q(\tau) = \frac{b_q + \bar{b}_{-q}}{\sqrt{2m\omega_q}}
