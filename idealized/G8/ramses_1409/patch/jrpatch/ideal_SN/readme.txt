This patch is written specifically for ideal supernovae simulations,
where a source emits radiation and winds for a mass-related lifetime before 
going SN..

Special considerations are:

TIMESTEP:
We override the large-by-default initial timestep (because the initial conditions
are 'static'). Also we make sure that SN blasts always happen at timestep 
intervals.

REFINEMENT:
We refine the region immediately around a SN source to the maximum
available level.
