#!/bin/bash

for i in {1..5}
do
   julia --threads 4 perturb_generator.jl $i 1 &
   julia --threads 4 perturb_generator.jl $i 2 &
   julia --threads 4 perturb_generator.jl $i 3 &
   julia --threads 4 perturb_generator.jl $i 4 &
   julia --threads 4 perturb_generator.jl $i 5 &
   julia --threads 4 perturb_generator.jl $i 6 &
   wait
done