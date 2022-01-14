% CALFEM - A Finite Element Toolbox.
% Version 3.6    2022-01-14
%
% Copyright (c)  Division of Structural Mechanics and
%                Divisuon of Solid Mechanics.
%                Lund University
%------------------------------------------------------------------------------
%
% Material functions.
%   hooke      - Form linear elastic constitutive matrix.
%   mises      - Compute stress and plastic strains for isotropic hardening von 
%                 Mises material.
%   dmises     - Form elasto-plastic contiuum matrix for isotropic hardening 
%                 von Mises material.
%
% Spring element functions.
%   spring1e   - Compute element matrix for spring element.
%   spring1s   - Compute spring force.
%
% Bar element functions.
%   bar1e      - Compute element matrix for 1D bar element.
%   bar1s      - Compute normal force in 1D bar element.
%   bar1we     - Compute element matrix for D1 bar element with elastic support.
%   bar1ws     - Compute normal force in 1D bar element with elastic support.
%   bar2e      - Compute element matrix for 2D bar element.
%   bar2s      - Compute normal force in 2D bar element.
%   bar2ge     - Compute element matrix for 2D geometrical nonlinear bar element.
%   bar2gs     - Compute normal force and axial force for 2D geometric nonlinear 
%                element.
%   bar3e      - Element matrix for 3D bar element.
%   bar3s      - Compute normal force in 3D bar element.
%
% Heat flow element functions.
%   flw2te     - Compute element matrix for triangular element.
%   flw2ts     - Compute temperature gradient and flux in triangular element.
%   flw2qe     - Compute element matrix for quadrilateral element.
%   flw2qs     - Compute temperature gradient and flux in quadrilateral element.
%   flw2i4e    - Compute element matrix for 4 node isoparametric element.
%   flw2i4s    - Compute temperature gradient and flux in 4 node isoparametric element.
%   flw2i8e    - Compute element matrix for 8 node 2D isoparametric element.
%   flw2i8s    - Compute temperature gradient and flux, 8 node 2D isoparametric element.
%   flw3i8s    - Compute element matrix for 8 node 3D isoparametric element.
%   flw3i8e    - Compute temperature gradient and flux in 8 node 3D isoparametric element.
%
% Solid element functions.
%   plante     - Compute element matrix for triangular element.
%   plants     - Compute stress and strain in triangular element.
%   plantf     - Compute internal element forces in triangular element.
%   planqe     - Compute element matrix for quadrilateral element.
%   planqs     - Compute stress and strain in quadrilateral element.
%   planre     - Compute element matrix for rectangular Melosh element.
%   planrs     - Compute stress and strain in rectangular Melosh element.
%   plantce    - Compute element matrix for rectangular Turner-Clough element.
%   plantcs    - Compute stress and strain in rectangular Turner-Clough element.
%   plani4e    - Compute element matrix for 4 node isoparametric element.
%   plani4s    - Compute stress and strain in 4 node isoparametric element.
%   plani4f    - Compute internal element forces in 4 node isoparametric element.
%   plani8e    - Compute element matrix for 8 node 2D isoparametric element.
%   plani8s    - Compute stress and strain in 8 node 2D isoparametric element.
%   plani8f    - Compute internal element forces, 8 node 2D isoparametric element.
%   soli8e     - Compute element matrix for 8 node 3D isoparametric element.
%   soli8s     - Compute stress and strain in 8 node 3D isoparametric element.
%   soli8f     - Compute internal element forces, 8 node 3D isoparametric element.
%
% Beam element functions.
%_  beam1e
%_  beam1we
%_  beam1we
%_  beam1ws
%   beam2e     - Compute element matrices for 2D beam element.
%   beam2s     - Compute section forces for 2D beam element.
%   beam2te    - Compute element matrices for 2D Timoshenko beam element.
%   beam2ts    - Compute section force, 2D Timoshenko beam element.
%   beam2we    - Compute element matrices for 2D beam element on elastic foundation.
%   beam2ws    - Compute section forces for 2D beam element on elastic foundation.
%   beam2ge    - Compute element matrix for geometrical nonlinear 2D beam element.
%   beam2gs    - Compute section forces for geometrical nonlinear 2D beam element.
%   beam2gxe   - Compute element matrix for geometrical nonlinear 2D exact beam element.
%   beam2gxs   - Compute section forces for geometrical nonlinear 2D exact beam element.
%   beam2de    - Compute element matrix for 2D beam element for dynamic analysis.
%   beam2ds    - Compute section forces for 2D beam element for dynamic analysis.
%   beam3e     - Compute element matrix for 3D beam element.
%   beam3s     - Compute section forces for 3D beam element.
%
% Plate element functions.
%   platre     - Compute element matrices for rectangular plate element.
%_  platrs     - Compute section forces for rectangular plate element.
%
% Static system functions.
%   assem      - Assemble element matrices.
%   coordxtr   - Extract element coordinates from a global coordinate matrix.
%   eigen      - Solve a generalized eigenvalue problem.
%   extract_ed - Extract values from a global vector.
%   insert     - Assemble element internal force vector.
%   red        - Reduce the size of a square matrix.
%   solveq     - Solve a system of equations.
%   statcon    - Perform static condensation.
%
% Dynamic system functions.
%   dyna2      - Solve a set of uncoupled second-order differential equations.
%   dyna2f     - Solve a set of uncoupled second-order differential equations 
%                in the frequency domain.
%   freqresp   - Compute frequency response.
%   gfunc      - Linear interpolation between equally spaced points.
%   ritz       - Approximative eigenvalues and eigenvectors by Lanczos method.
%   spectra    - Compute seismic response spectra.
%   step1      - Step-by-step integration in first-order systems.
%   step2      - Step-by-step integration in second-order systems.
%   sweep      - Frequency response function.
%
% Graphics functions.
%   dispbeam2  - Draw displacements for beam element.
%   eldraw2    - Draw undeformed finite element mesh.
%   eldisp2    - Draw deformed finite element mesh.
%   elflux2    - Plot flux vectors.
%   eliso2     - Draw isolines for nodal quantities.
%   elprinc2   - Plot principal stresses.
%   scalfact2  - Determine scale factor.
%   scalgraph2 - Draw graphic scale.
%   secforce2  - Draw section force diagram for bar or beam element.
%