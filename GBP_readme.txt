Requirements for rotating reference frame and rotationally periodic BCs:

1) Translational periodic BCs no longer allowed.
2) Parameter 'ldg-beta' in [solver-interfaces] must be set to 0.5 if running Navier-Stokes.
3) Parameter 'omg' in [solver-constants] must be defined giving the angular frequency (radians per time unit).
   Setting 'omg' = 0.0 will turn off the rotating frame terms. 
   Angular frequency vector is aligned with +z axis (omg > 0 -> counter-clockwise rotation viewed from above).
4) Flow/geometry must be aligned such that axis extending out of the top of the rotor is +z.
   Rotation acts on the x/y plane (u/v-velocity components). z (w-velocity) is unaffected.