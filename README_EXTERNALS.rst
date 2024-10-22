Obtaining optional external libraries that can be used with mizuRoute
=====================================================================

mizuRoute is obtained via github. And optionally there are some external libraries
that can be used with it. These libraries are obtained by git-fleximod utility


To obtain the mizuRoute code and the optional externals you need to do the following:

#. Clone the repository. ::

      git clone https://github.com/ESCOMP/mizuRoute.git my_mizuRoute_sandbox

   This will create a directory ``my_mizuRoute_sandbox/`` in your current working directory.

#. Run **git-fleximod**. ::

      ./bin/git-fleximod -g .gitmodules update 

   **./bin/git-fleximod** is a python wrapper script that calls actual git-fleximod located under .lib/git-fleximod. 
   The git-fleximod ulitity downloads source codes of external libraries specified in **.gitmodules**. mizuRoute now uses 
   only one external library (parallel-io).

   NOTE: 

   This second step is only required if you are going to use
   one of the external libraries in the mizuRoute build (currently pio).

   If this is the first time to download external libraries, you may need to create a directory under libraries 
   before run **git-fleximod**. e.g.,::

     mkdir -p libraries/parallelio

At this point you have a working version of the full mizuRoute and all optional libraries.

More details on checkout_externals
----------------------------------

The file **.gitmodules** in your top-level mizuRoute directory tells
**git-fleximod** which tag/branch of external libraries
that should be brought in to generate your sandbox.

NOTE: git-fleximod will always attempt to make the working copy 
exactly match the externals description. If
you manually modify an external without updating .gitmodules, e.g. switch
to a different tag, then rerunning checkout_externals will switch you
back to the external described in .gitmodules. See below
documentation `Customizing your mizuRoute sandbox`_ for more details.

**You need to rerun ./bin/git-fleximod whenever .gitmodules has
changed** (unless you have already manually updated the relevant
external(s) to have the correct branch/tag checked out). Common times
when this is needed are:

* After checking out a new mizuRoute branch/tag

* After merging some other mizuRoute branch/tag into your currently
  checked-out branch

**checkout_externals** must be run from the root of the source
tree. For example, if you cloned mizuRoute with::

  git clone https://github.com/ESCOMP/mizuRoute.git my_mizuRoute_sandbox

then you must run **checkout_externals** from
``/path/to/my_mizuRoute_sandbox``.

To see more details of **git-fleximod**, issue ::

  ./bin/git-fleximod --help

Customizing your mizuRoute sandbox
==================================

There are several use cases to consider when you want to customize or modify your mizuRoute sandbox.

Switching to a different mizuRoute branch or tag
-------------------------------------------

If you have already checked out a branch or tag and **HAVE NOT MADE ANY
MODIFICATIONS** it is simple to change your sandbox. Say that you
checked out cesm-coupling.n02_v2.1.2 but really wanted to have cesm-coupling.n02_v2.1.3;
you would simply do the following::

  git checkout cesm-coupling.n02_v2.1.3
  ./bin/git-fleximod -g .gitmodules update

You should **not** use this method if you have made any source code
changes, or if you have any ongoing mizuRoute simulations that you want
to saze that were created from this sandbox. In these cases, it is often 
easiest to do a second **git clone**.

Pointing to a different version of an external
----------------------------------------------

Each entry in **.gitmodules** has the following form
::
  [submodule "parallelio"]
  path = libraries/parallelio
  url = https://github.com/NCAR/ParallelIO
  fxtag = pio2_6_2
  fxrequired = ToplevelRequired
  # Standard Fork to compare to with "git fleximod test" to ensure personal forks aren't committed
  fxDONOTUSEurl = https://github.com/NCAR/ParallelIO

**fxtag** entry specifies either a tag or a branch. To point to a new tag:

#. Modify **fxtag** in **.gitmodules** (e.g., changing
   ``pio2_6_2`` to ``pio2_6_3`` above)

#. Update a tag of the library::

     ./bin/git-fleximod -g .gitmodules update

Keep in mind that randomly changing individual externals from a tag may result
in an invalid model (won't compile or won't run).
However, since the external components being used are robust changing
them is unlikely to have an effect outside of possibly interfering with
the ability to build with the externals.

Committing your change to .gitmodules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After making this change, it may be a good idea to commit the change in your
local mizuRoute git repository. First create a branch in your local
repository, then commit it. For example::

  git checkout -b my_mizuRoutebranch
  git add .gitmodules
  git commit -m "Update PIO to pio2_6_3"


