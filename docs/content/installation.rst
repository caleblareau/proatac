.. _installation:

Installation
============

BEDOPS is available to users as :ref:`pre-built binaries <installation_via_packages>` and :ref:`source code <installation_via_source_code>`.

.. _installation_via_packages:

======================
Via pre-built packages
======================

Pre-built binaries offer the easiest and fastest installation option for users of BEDOPS. At this time, we offer binaries for 32- and 64-bit versions of Linux and OS X (Intel) platforms.

-----
Linux
-----

1. Download the current 32- or 64-bit package for Linux from `Github BEDOPS Releases <https://github.com/bedops/bedops/releases>`_.
2. Extract the package to a location of your choice. 
   In the case of 32-bit Linux: ::

       $ tar jxvf bedops_linux_i386-vx.y.z.tar.bz2

   In the case of 64-bit Linux: ::

       $ tar jxvf bedops_linux_x86_64-vx.y.z.tar.bz2

   Replace ``x``, ``y`` and ``z`` with the version number of BEDOPS you have downloaded.
3. Copy the extracted binaries to a location of your choice which is in your environment's ``PATH``, *e.g.* ``/usr/local/bin``: ::

       $ cp bin/* /usr/local/bin

   Change this destination folder, as needed.

--------
Mac OS X
--------

1. Download the current Mac OS X package for BEDOPS from `Github BEDOPS Releases <https://github.com/bedops/bedops/releases>`_.
2. Locate the installer package (usually located in ``~/Downloads`` |--| this will depend on your web browser configuration):

   .. image:: ../assets/installation/bedops_macosx_installer_icon.png

3. Double-click to open the installer package. It will look something like this:

   .. image:: ../assets/installation/bedops_macosx_installer_screen_v2.png
      :width: 99%

4. Follow the instructions to install BEDOPS and library dependencies to your Mac. (If you are upgrading from a previous version, components will be overwritten or removed, as needed.)

.. _installation_via_source_code:

===============
Via source code
===============

.. _installation_via_source_code_on_linux:

-----
Linux
-----

Compilation of BEDOPS on Linux requires GCC 4.8.2 (both ``gcc`` and ``g++`` and related components) or greater, which includes support for `C++11 <http://en.wikipedia.org/wiki/C%2B%2B11>`_ features required by core BEDOPS tools. Other tools may be required as described in the installation documentation that follows.

.. |--| unicode:: U+2013   .. en dash
.. |---| unicode:: U+2014  .. em dash, trimming surrounding whitespace
   :trim:
