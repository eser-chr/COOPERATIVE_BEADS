# Beads_on_line
Monte Carlo simulation, TASEP with reaction terms


The code language is C++ and the code was created for my master's thesis.
The link for the script is: https://studenttheses.uu.nl/bitstream/handle/20.500.12932/39774/Thesis_Chris.pdf?sequence=1&isAllowed=y



Initial setup:

I created a 1D discrete road. Either randomly or based on a given distribution i placed particles on the axis.



Dynamics:


For each step the programm will shuffle the order of particles to perform their move. This is done to avoid correlations that would appear if we had a queue of particles, that each of them would perform the following.


Each bead can be found into two states. The one we can call it bound and the other unbound.


BOUND STATE: The particles have a probability to step forward or to unbind. Using a random number we can simulate this mechanism. 
If it tries to move but the position is obtained by another particle the step is cancelled.
If it tries to unbind, it will unbind.
If it tries to do nothing, nothing will change.


UNBOUND STATE: The particles can only bind to the road. We simulate this as previously with a random number.



Output:

The script writes in a file all the positions. We can analyze them later on with python.
A simple example for how to analyze is the file plotter.py.
