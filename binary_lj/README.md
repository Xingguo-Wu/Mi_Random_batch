# Kob–Andersen binary LJ

`in.binary_lj` is the Mi-RBL input for the 80:20 mixture. Its defaults give the paper case: 8,788 particles, density 1.20, temperature ramp from 1.00 to 0.60, `P=10`, `mi_alpha=0.31`, 150,000 equilibration steps, and 400,000 production steps.

The eight seeds are `12345`, `22345`, `32345`, `42345`, `52345`, `62345`, `72345`, and `82345`. `run_example.sh` provides a short check.

For comparison runs, copy the input and change the method commands:

| Method | Changes |
| --- | --- |
| RBL | Set `-var mi_alpha 0`. |
| CL | Use `pair_style lj/cut ${cutmax}`; remove `fix mirbl_state` and the gain dump. |
| VR-RBL | Use `pair_style lj/cut/vrrbl ${cutmax} ${rc} ${p} ${seed}`; replace both Langevin fixes with `vr/langevin`; remove `fix mirbl_state` and the gain dump. |

The log contains kinetic energy; `rdf.binary_lj_ka.dat`, `traj_unwrapped.lammpstrj`, and `pro_mirbl_state.dump` contain RDFs, coordinates, and gains. The paper uses production steps 150,000–350,000 for temperature analysis and 200,000–400,000 for RDFs.
