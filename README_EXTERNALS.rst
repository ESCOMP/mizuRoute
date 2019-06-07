Obtaining optional external libraries that can be used with mizuRoute
=====================================================================

mizuRoute is obtained via github. And optionally there are some external libraries
that can be used with it. These libraries are obtained by r


To obtain the mizuRoute code and the optional externals you need to do the following:

#. Clone the repository. ::

      git clone https://github.com/NCAR/mizuRoute.git my_mizuRoute_sandbox

   This will create a directory ``my_mizuRoute_sandbox/`` in your current working directory.

#. Run the script **manage_externals/checkout_externals**. ::

      ./manage_externals/checkout_externals

   The **checkout_externals** script is a package manager that will
   populate the mizuRoute directory with the CIME infrastructure code that can
   then be used. 

   NOTE: This second step is only required if you are going to use
   one of the external libraries in the mizuRoute build (currently either mpi-serial, or pio).

At this point you have a working version of the full mizuRoute and all optional libraries.

More details on checkout_externals
----------------------------------

The file **Externals.cfg** in your top-level mizuRoute directory tells
**checkout_externals** which tag/branch of cime (and possibly other externals)
that should be brought in to generate your sandbox.

NOTE: checkout_externals will always attempt to make the working copy 
exactly match the externals description. If
you manually modify an external without updating Externals.cfg, e.g. switch
to a different tag, then rerunning checkout_externals will switch you
back to the external described in Externals.cfg. See below
documentation `Customizing your mizuRoute sandbox`_ for more details.

**You need to rerun checkout_externals whenever Externals.cfg has
changed** (unless you have already manually updated the relevant
external(s) to have the correct branch/tag checked out). Common times
when this is needed are:

* After checking out a new mizuRoute branch/tag

* After merging some other mizuRoute branch/tag into your currently
  checked-out branch

**checkout_externals** must be run from the root of the source
tree. For example, if you cloned mizuRoute with::

  git clone https://github.com/NCAR/mizuRoute.git my_mizuRoute_sandbox

then you must run **checkout_externals** from
``/path/to/my_mizuRoute_sandbox``.

To see more details of **checkout_externals**, issue ::

  ./manage_externals/checkout_externals --help

Customizing your mizuRoute sandbox
==================================

There are several use cases to consider when you want to customize or modify your mizuRoute sandbox.

Switching to a different mizuRoute branch or tag
-------------------------------------------

If you have already checked out a branch or tag and **HAVE NOT MADE ANY
MODIFICATIONS** it is simple to change your sandbox. Say that you
checked out v1.0.1 but really wanted to have v1.0.2;
you would simply do the following::

  git checkout v1.0.2
  ./manage_externals/checkout_externals

You should **not** use this method if you have made any source code
changes, or if you have any ongoing mizuRoute simulations that you want
to saze that were created from this sandbox. In these cases, it is often 
easiest to do a second **git clone**.

Pointing to a different version of an external
----------------------------------------------

Each entry in **Externals.cfg** has the following form
::

  [cime]
  local_path = cime
  protocol = git
  repo_url = https://github.com/ESMCI/cime
  tag = cime5.8.2
  required = True

Each entry specifies either a tag or a branch. To point to a new tag:

#. Modify the relevant entry/entries in **Externals.cfg** (e.g., changing
   ``cime5.8.0`` to ``cime5.8.2`` above)

#. Checkout the new component(s)::

     ./manage_externals/checkout_externals

Keep in mind that randomly changing individual externals from a tag may result
in an invalid model (won't compile or won't run).
However, since the external components being used are robust changing
them is unlikely to have an effect outside of possibly interfering with
the ability to build with the externals.

Committing your change to Externals.cfg
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After making this change, it's a good idea to commit the change in your
local CTSM git repository. First create a branch in your local
repository, then commit it. For example::

  git checkout -b my_mizuRoutebranch
  git add Externals.cfg
  git commit -m "Update CIME to cime5.8.2"


