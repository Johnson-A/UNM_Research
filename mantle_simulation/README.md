# Mantle Simulation
Simulation for mantle convection with fluid flow and heat advection written in
FENICS. There are two main components for the code, the simulation itself and
the post processing.

## Simulation
The code consists of the following files

* `__init__.py`

This is for now simply an empty file which marks this directory as one
containing a package.
http://stackoverflow.com/questions/448271/what-is-init-py-for

* `constants.py`

Contains some of the simulation constants like mesh width and height, and mesh
fineness. Additional physical parameters could be moved here as well to simplify
simulation.py.

* `LAB.py`

Definition file for the lithosphere-asthenosphere boundary. The profile is
generated using two symmetric and reflected hyperbolic tangent functions. The
variable parameters for the profile include the keel_width which defines the
mirrored offset from the center of the mesh about which both tanh are
calculated. The LAB_height defines the background LAB height. The steepness of
the keel sides are determined by the scale. The keel_height is the depth of the
keel.

* `simulation.py`

Here is where the simulation is actually implemented. It is best to read through
the code and the corresponding comments, but a high level overview is as follows.

Setup multi-threaded organization through the FENICS MPI interface. The
main_proc decorator is used to ensure only the master thread with RANK == 0
executes certain statements. Some synchronization is done when creating the
output directory.

Set form_compiler options for performance.

Calculate physical constants in dimensionless quantities.

Define PeriodicBoundary conditions, Initial conditions, and the mesh and
associated function spaces. Define forms and time evolution scheme. Solve
and output files in main time evolution loop.


### Key Details
* XDMF File Format

We use the xdmf file format for simulation output because of its support for
concurrent output, viewing, and efficiency. Some references
http://www.xdmf.org/index.php/XDMF_Model_and_Format
https://support.hdfgroup.org/HDF5/

*Note* that when loading these files in a program like paraview, select the xdmf
file and not the hdf5 file which contains the actual data. There may be
additional load options for paraview, like `xdmf3`, but you should select the
generic `xdmf` or else paraview will crash.

* Creating and using new files

The DefaultDictByKey class extends from the builtin defaultdict collection in
python. This is an ordinary dictionary in python except when a key is accessed,
the entry is created with a function that takes the index as a parameter if it
does not already exists. This function in the code is create_xdmf. This means
that you don't have to explicitly setup the file and remember a variable name.
To create a new file, simply index with a new file name and you will be able to
write to it immediately.

*Note* In order to play the simulation back in paraview, all writes to a
single file must come from the *same* function in FENICS. If an intermediary
function is created each time, e.g. because one is not created and then
reassigned to, then they may be interpreted as separate data series.

* Notes

There is a file notes.org which keeps track of potential code changes and other
information.

### Running
```mpirun -n $num_threads python simulation.py $output_directory``` This will
run the simulation in parallel using $num_threads. To avoid completely
overwhelming the system, prefix the above command with `nice` or reduce
$num_threads. The optimal value for lithos would be num_threads = 8. When
running the code, consider using tmux to create a local user session and then
executing the above command. If the network is disconnected or you need to close
ssh, you can reattach to that session. A good guide is
https://danielmiessler.com/study/tmux/#gs.zMou=tA. The output directory can be
specified with the first commandline argument to simulation.py. When no argument
is passed, the default is 'run'. When setting up the base directory, if the
desired directory already exists, that directory will be backed up with its
current name and timestamp. A directory called code_copy will be created in the
output directory and contains the repository sha hash for the current git HEAD
and a list of all changes and their diffs. This makes it easy to verify what
code generated the results.


## Post-Processing
All post processing scripts are inside post_process/. The original 3-dimensional
code which calculated things like the dynamic pressure gradients at the surface
are in "StreamlineTplot.m" while the modified code for the 2d advection
calculation is in "StreamlineTplot_2d.m". The altmany-export-fig directory
contains the popular `export_fig` code available at
https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig and is
used to actually render the figures to a png file. *You will need to add the
whole repository to your matlab path because some code is reused between
different projects, like map_columns*. You will need to set the root_dir for the
output directory you wish to visualize. Figures are generated in groups, in
post_process/figures/ according to what they depict.

## Schumann
The original implementation of the Schumann model is a matlab program inside schumann.m

## Tests
These are miscellaneous scripts inside the tests/.


