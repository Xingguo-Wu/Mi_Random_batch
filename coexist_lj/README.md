# Gas–liquid coexistence LJ

`in.coexist_lj` is the Mi-RBL input. Its defaults give 7,000 particles at density 0.2 and temperature 0.9, with `P=5`, `alpha=0.59`, 200,000 equilibration steps, and 300,000 production steps.

The eight seeds are `12345`, `22345`, `32345`, `42345`, `52345`, `62345`, `72345`, and `82345`. The paper uses `alpha=0.58,0.59,0.60` for `P=5`; `0.36,0.37,0.38` for `P=10`; and `0.26,0.27,0.28` for `P=15`. Set these with `-var shell_samples P -var alpha ALPHA`. `run_example.sh` provides a short check.

For comparison runs, copy the input and change the method commands:

| Method | Changes |
| --- | --- |
| RBL | Set `-var alpha 0`. |
| CL | Set `newton on`, use `pair_style lj/cut ${cut_lj}`, and remove `fix mirbl_state` and the gain dump. |
| VR-RBL | Use `pair_style lj/cut/vrrbl ${cut_lj} ${core_size} ${shell_samples} ${seed}` and replace the Langevin fix with `vr/langevin`; remove `fix mirbl_state` and the gain dump. |
| Shared-gain ablation | Use `pair_style lj/cut/mirbl/ablation/global`; remove the gain dump. |

`traj_unwrapped.lammpstrj` contains coordinates and `pro_mirbl_state.dump` contains gains. For the paper MSD, use origins 400,000–490,000 every 10,000 steps, lags 0–8,000 steps, remove center-of-mass drift, and average over the eight seeds.
