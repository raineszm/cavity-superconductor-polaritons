.. Bardasis Schrieffer Polaritons documentation master file, created by
   sphinx-quickstart on Wed May 16 11:45:15 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Bardasis Schrieffer Polaritons's documentation!
==========================================================

.. _units:

Choice of units
---------------

Calculations performed using `atomic units <https://en.wikipedia.org/wiki/Atomic_units>`_ with the energy unit as `meV` instead of `Hartrees`.


.. _xi-approx:

xi approximation
-----------------

Several integrals are performed within the :math:`\xi` approximation. Specifically given an integral

.. math::

    \int \frac{d\mathbf{k}}{(2\pi)^2} f(\xi(\mathbf k), \mathbf k)

we rewrite it as

.. math::

    \nu \int_0^\infty d\xi \int_0^{2\pi}\frac{d\theta}{2\pi} \left[f(\xi, \mathbf{k}_f) + f(-\xi, \mathbf{k}_f)\right]

For the case of a function which depends only on the quasparticle energy :math:`\lambda_k = \sqrt{\xi^2 + \Delta^2}` this can also be written

.. math::

    2 \nu \int_\Delta^\infty d\lambda \frac{\lambda}{\sqrt{\lambda^2 - \Delta^2}}\int_0^{2\pi}\frac{d\theta}{2\pi} f(\lambda, \mathbf{k}_f)]

Physics
=======

.. toctree::
    :caption: Physics concepts

    modes


Code
====



C++
----

.. toctree::
   :maxdepth: 2
   :caption: C++ modules:

   system
   state
   bs
   cavity
   coupling
   polariton

Python
------

.. toctree::
    :caption: Python modules
    :glob:

    pymodules/*



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
