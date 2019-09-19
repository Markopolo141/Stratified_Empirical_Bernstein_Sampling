# Stratified_Empirical_Bernstein_Sampling
Source for the paper "Stratified finite Empirical Bernstein Sampling"
compile all data used in the paper by executing run.sh script.
required is a version of Python 2.7, and gnu C++ compiler under posix compliant system.
optionally install the python library 'tqdm' for some pretty progress bars.

Please Note: on the 28th of August, two implimentations bugs were found,
The first didnt affect the results but resulted in periodic stalling of the computation,
The second bug was about the implementation of Neyman sampling (not in the Shapley value case)
and fixing it made our results seem slightly less impressive in that context.
We plan to bring it up with the editor.
