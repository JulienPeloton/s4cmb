This folder contains example scripts to run basic and advanced App.

test
===============
This folder contains a simple App to test the code. It is also used
by Travis to check the code after each push. You can test the code by
running:

::

    python examples/test/simple_app.py -inifile examples/inifiles/simple_parameters.py -tag 'test'

perfect
===============
This folder contains an App to test the library on real condition (large fraction)
of sky, large number of detector, and so on. but without systematic effects.
You can launch it on cluster using the provided launcher (made for NERSC) in
the same directory. You might need to update the paths in the launcher to
match your tree.

contamination
===============
