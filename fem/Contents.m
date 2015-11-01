% CALFEM - A Finite Element Toolbox.
% Version 3.3  99-03-01
%
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%------------------------------------------------------------------------------
%
% Matrix functions.
%   red        - Reduce the size of a square matrix.
%
% Static system functions.
%   assem      - Assemble element matrices.
%   coordxtr   - Extract element coordinates from a global coordinate matrix.
%   eigen      - Solve a generalized eigenvalue problem.
%   extract    - Extract values from a global vector.
%   insert     - Assemble element internal force vector.
%   solveq     - Solve a system of equations.
%   statcon    - Perform static condensation.
%
% Dynamic system functions.
%   dyna2      - Solve set of uncoupled second-order diff. eqs.
%   dyna2f     - Solve set of uncoupled second-order diff. eqs. in freq. domain.
%   freqresp   - Compute frequency response.
%   gfunc      - Linear interpolation between equally spaced points.
%   ritz       - Approximative eigenvalues and eigenvectors by Lanczos method.
%   spectra    - Compute seismic response spectra.
%   step1      - Step-by-step integration in first-order systems.
%   step2      - Step-by-step integration in second-order systems.
%   sweep      - Frequency response function.
%
% Spring element functions.
%   spring1e   - Element matrix for spring element.
%   spring1s   - Compute spring force.
%
% Bar element functions.
%   bar2e      - Element matrix for 2D bar element.
%   bar2g      - Element matrix for 2D geometrical nonlinear bar element.
%   bar2s      - Compute normal force in 2D bar element.
%   bar3e      - Element matrix for 3D bar element.
%   bar3s      - Compute normal force in 3D bar element.
%
% Beam element functions.
%   beam2e     - Element matrix, 2D beam element.
%   beam2s     - Section force, 2D beam element.
%   beam2t     - Element matrix, 2D Timoshenko beam element.
%   beam2ts    - Section force, 2D Timoshenko beam element.
%   beam2w     - Element matrix, 2D beam element on elastic foundation.
%   beam2ws    - Section force, 2D beam element on elastic foundation.
%   beam2g     - Element matrix, geometrical nonlinear 2D beam element.
%   beam2gs    - Section force, geometrical nonlinear 2D beam element.
%   beam2d     - Element matrix, 2D beam element for dynamic analysis.
%   beam2ds    - Section force, 2D beam element for dynamic analysis.
%   beam3e     - Element matrix, 3D beam element.
%   beam3s     - Section force, 3D beam element.
%
% Heat flow element functions.
%   flw2te     - Element matrix, triangular element.
%   flw2ts     - Temperature gradient and flux, triangular element.
%   flw2qe     - Element matrix, quadrilateral element.
%   flw2qs     - Temperature gradient and flux, quadrilateral element.
%   flw2i4e    - Element matrix, 4 node isoparametric element.
%   flw2i4s    - Temperature gradient and flux, 4 node isoparametric element.
%   flw2i8e    - Element matrix, 8 node 2D isoparametric element.
%   flw2i8s    - Temperature gradient and flux, 8 node 2D isoparametric element.
%   flw3i8s    - Element matrix, 8 node 3D isoparametric element.
%   flw3i8e    - Temperature gradient and flux, 8 node 3D isoparametric element.
%
% Solid element functions.
%   plante     - Element matrix, triangular element.
%   plants     - Stress and strain, triangular element.
%   plantf     - Internal element forces, triangular element.
%   planqe     - Element matrix, quadrilateral element.
%   planqs     - Stress and strain, quadrilateral element.
%   planre     - Element matrix, rectangular Melosh element.
%   planrs     - Stress and strain, rectangular Melosh element.
%   plantce    - Element matrix, rectangular Turner-Clough element.
%   plantcs    - Stress and strain, rectangular Turner-Clough element.
%   plani4e    - Element matrix, 4 node isoparametric element.
%   plani4s    - Stress and strain, 4 node isoparametric element.
%   plani4f    - Internal element forces, 4 node isoparametric element.
%   plani8e    - Element matrix, 8 node 2D isoparametric element.
%   plani8s    - Stress and strain, 8 node 2D isoparametric element.
%   plani8f    - Internal element forces, 8 node 2D isoparametric element.
%   soli8e     - Element matrix, 8 node 3D isoparametric element.
%   soli8s     - Stress and strain, 8 node 3D isoparametric element.
%   soli8f     - Internal element forces, 8 node 3D isoparametric element.
%
% Plate element functions.
%   platre     - Element matrix, rectangular plate element.
%   platrs     - Section forces, rectangular plate element.
%
% Material functions.
%   hooke      - Linear elastic constitutive matrix.
%   mises      - Stress and plastic strain, isotropic hardening von Mises mtrl.
%   dmises     - Elasto-plastic matrix, isotropic hardening von Mises mtrl.
%
% Graphics functions.
%   eldraw2    - Draw undeformed finite element mesh.
%   eldisp2    - Draw deformed finite element mesh.
%   eldia2     - Draw section force diagram.
%   elflux2    - Plot flux vectors.
%   eliso2     - Draw isolines for nodal quantities.
%   elprinc2   - Plot principal stresses.
