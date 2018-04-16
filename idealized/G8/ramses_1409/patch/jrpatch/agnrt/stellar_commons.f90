module stellar_commons
   integer, parameter::nbint=100
   integer, parameter::nbinz=5
   real(KIND=8),dimension(1:nbint)        :: t_wind
   real(KIND=8),dimension(1:nbinz)        :: Z_wind
   real(KIND=8),dimension(1:nbint,1:nbinz):: cMwind    !cumulative mass
   real(KIND=8),dimension(1:nbint,1:nbinz):: cEwind    !cumulative energy
   real(KIND=8),dimension(1:nbint,1:nbinz):: cMZwind   !cumulative mass*metal
   real(KIND=8),dimension(1:nbint,1:nbinz):: cMOxwind  !cumulative mass*O
   real(KIND=8),dimension(1:nbint,1:nbinz):: cMFewind  !cumulative mass*Fe
end module stellar_commons
