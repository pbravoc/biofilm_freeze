#!/bin/bash

for i in {1..8}
do
   julia --threads 4 u0_generator.jl $i 1 &
   julia --threads 4 u0_generator.jl $i 2 &
   julia --threads 4 u0_generator.jl $i 3 &
   julia --threads 4 u0_generator.jl $i 4 &
   julia --threads 4 u0_generator.jl $i 5 &
   wait
done