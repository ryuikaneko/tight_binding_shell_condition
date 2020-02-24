# check shell conditions of tight binding models

* Check open or closed shell for each system (closed shell recommended for QMC)

* Lattices
  * 1D
    * Chain: E(k) = -2cos(k)
  * 2D
    * Square: E(kx,ky) = -2(cos(kx)+cos(ky))
    * Triangular: E(kx,ky) = -2(cos(kx)+cos(ky)+cos(kx+ky))
    * Honeycomb: E_{1,2}(kx,ky) = \pm |1 + exp(-i kx) + exp(-i ky)|
    * Kagome: E_{1,2,3}(kx,ky) = 2, -1 \pm \sqrt(1 + 8 cos(kx/2) cos(ky/2) cos((kx-ky)/2))

* Boundary
  * 1D
    * Periodic (P)
    * Antiperiodic (AP)
  * 2D
    * P-P
    * P-AP
    * AP-P
    * AP-AP

* In data files,
  * "Energy density" = E/Nsite/Nspin
  * "Number of electrons" = number of up spins = number of down spins

* Results (DOS)
  * Chain
![DOS chain](https://raw.githubusercontent.com/ryuikaneko/tight_binding_shell_condition/master/1d_chain/filling_1over2_BC_AP/fig_1d_chain_dos.png "DOS chain")

  * Square
![DOS square](https://raw.githubusercontent.com/ryuikaneko/tight_binding_shell_condition/master/2d_square/filling_1over2_BC_P_AP/fig_2d_square_dos.png "DOS square")

  * Triangular
![DOS triangular](https://raw.githubusercontent.com/ryuikaneko/tight_binding_shell_condition/master/2d_triangular/filling_1over2_BC_P_AP/fig_2d_triangular_dos.png "DOS triangular")

  * Honeycomb
![DOS honeycomb](https://raw.githubusercontent.com/ryuikaneko/tight_binding_shell_condition/master/2d_honeycomb/filling_1over2_BC_P_AP/fig_2d_honeycomb_dos.png "DOS honeycomb")

  * Kagome
![DOS kagome](https://raw.githubusercontent.com/ryuikaneko/tight_binding_shell_condition/master/2d_kagome/filling_1over2_BC_P_AP/fig_2d_kagome_dos.png "DOS kagome")
