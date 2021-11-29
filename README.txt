File structure

main.m : main file, all other functions are called from this, not parallelized
main_fast.m: same as main.m but it is parallelized and optimized to run faster

Input.m: Input file, includes all parameters required for the simulation

transmission: function that calculates transmission
transmission_alt.m: parallelized and optimized function for transmission



