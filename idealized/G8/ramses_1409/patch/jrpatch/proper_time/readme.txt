Patch to keep track of the proper age of particles. 

A new property, tpp, has been added to particles, signifying their proper birthtime.

Using this property, the subroutine 
 	getAgeGyr(tpp(ind_part), age)    (in update_time.f90)
can be used to give the proper age of a particle.

Joki Rosdahl, Sept 2010.