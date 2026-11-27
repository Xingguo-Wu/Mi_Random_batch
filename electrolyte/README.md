# Primitive electrolyte

`in.electrolyte` is the Mi-RBE input for a 3:1 electrolyte. Its defaults give 800 cations and 2,400 anions in a box of side 25.6, with `P=10`, `mi_alpha=0.4`, 150,000 equilibration steps, and 200,000 production steps.

The eight seeds are `11031`–`11038`. Set `P` and `alpha` with `-var p P -var mi_alpha ALPHA`. For the charge-density figure, use `alpha=0.4,0.5,0.6` at `P=5,10` and `alpha=0.1,0.2,0.3` at `P=20`. The kinetic and gain panels use `P=5,10` with `alpha=0.4,0.5,0.6`; the timing panel uses `P=20, alpha=0.2`. For timing runs, keep the ion ratio 1:3 and scale the box side as `25.6*(N/3200)^(1/3)`. `run_example.sh` provides a short check.

For comparison runs, copy the input and change the method commands:

| Method | Changes |
| --- | --- |
| RBE | Set `-var mi_alpha 0`. |
| PPPM | Use `kspace_style pppm 1.0e-4`; remove `fix mirbe_state`, `fix mirbe_gain`, and the gain dump. |

The log contains kinetic energy; `traj.electrolyte_bulk.lammpstrj` and `mirbe_state.dump` contain coordinates and gains. Processed figure data are in `../paper/data/`.
