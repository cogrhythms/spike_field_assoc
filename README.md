spike_field_assoc
=================

A fork of Kyle LePage's spike field association code.

The parent of this source is the code provided on [Kyle LePage's website](http://math.bu.edu/people/lepage/code.html), dated 20 Nov 2012.


Setup instructions
------------------

Run get_cases_param_est_with_exp_many_f.m

The .m files:
glm_sfa_*_chk_phase.m  produce 2-page postscript files with diagnostics.

The first page is a plot of the filter transfer function and the second plot (right most column) has a plot of instantaneous frequency departure.

These plots should be "green everywhere" except at the trial boundaries.  If the n_trim parameter is mis-specified you will see a big departure from
this behaviour; or if there is a problem with "leakage from adjacent frequencies".

