# check shell conditions of tight binding models

* Check open or closed shell for each system (closed shell recommended for QMC)
* Lattices
  * 1D
    * Chain: E(k) = -2cos(k)
  * 2D
    * Square: E(kx,ky) = -2(cos(kx)+cos(ky))
    * Triangular: E(kx,ky) = -2(cos(kx)+cos(ky)+cos(kx+ky))
    * Honeycomb: E_{1,2}(kx,ky) = \pm |1+exp(-i*kx)+exp(-i*ky)|
    * Kagome: E_{1,2,3}(kx,ky) = 2, -1 \pm \sqrt(1+8*cos(kx/2)*cos(ky/2)*cos((kx-ky)/2))
* Boundary
  * 1D
    * Periodic (P)
    * Antiperiodic (AP)
  * 2D
    * P-P
    * P-AP
    * AP-P
    * AP-AP
  
