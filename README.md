# FluidSimulation
In this project I implemented a fluid simulator using FLIP method in C++ using OpenVDB 9.

## Implemented Features
* FLIP/PIC method
* 3D simulation
* OpenVDB integration
* Finite Volumes method to better represent fluid-solid and fluid-air interfaces
* Strong two-way coupling
* Density Projection to preserve the fluid volume. This is my implementation of [this paper](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8869736) by Tassilo Kugelstadt, Andreas Longva, Nils Thuerey, and Jan Bender

## Demonstration
The videos below are using OpenVDB's `vdb_view` as a visualisation tool. This is not a render just a simple visualisation of the zero level set.

- Simple drop of fluid cube

https://user-images.githubusercontent.com/32745362/195076875-d98e46e0-6056-4fd8-a3c2-174f9ed9f12b.mov

- Side drop of fluid cube

https://user-images.githubusercontent.com/32745362/195076624-5b57fb02-96fc-4317-8ad2-05072a138c03.mov

- Ball of half the water's density floating

https://user-images.githubusercontent.com/32745362/195076931-31b83383-8198-4468-b274-fef41404b7e5.mov

- Ball of half the water's density pushed sideways

https://user-images.githubusercontent.com/32745362/195076983-ce7dfc69-d9ae-4bc5-8657-08c11952df4f.mov

