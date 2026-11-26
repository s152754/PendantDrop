# Task: reorganize uncertainty functions and write an overview document

Scope: the project scripts `Camera_intrinsics.jl`, `Camera_extrinsics.jl` and `uncertainty_camtodrop.jl`,
and the source file `camera_uncertainties.jl`. Do not modify `camera.jl` or any other file.

General rules:
- Keep the `abspath(PROGRAM_FILE) == @__FILE__ || isinteractive()` guard in all three scripts, so I can
  still run them with the VS Code run button.
- Run and verify everything in the project environment (`--project=.`).
- If something is unclear, ask rather than guess.

## Phase 1 — analysis and move plan (change nothing yet)

1. **Classify functions.** For every function defined in the three scripts, decide whether it is
   uncertainty-related (computes, samples or propagates parameter uncertainty: bootstrap procedures,
   covariance estimation, Jacobians for error propagation, Monte Carlo sampling of camera parameters)
   or not (camera calibration, ArUco detection, pose estimation, image loading, plotting, I/O).
   Uncertainty-related functions will move to `camera_uncertainties.jl`; all others stay where they are.
   List ambiguous cases separately (e.g. a bootstrap function that calls a pose solver: the bootstrap
   moves, the pose solver stays) and explain your proposed choice.
2. **Check dependencies.** For each function to move, report whether it depends on variables, constants
   or functions that exist only in its script. If moving it would require more than updating calls and
   imports, flag it instead of proposing a workaround.
3. **Module structure.** `camera_uncertainties.jl` is accessed as `Model.Camera.IntrinsicUncertainty`.
   Explain how moved functions will be made accessible from the scripts (exports or module prefixes).
   The module name no longer matches its contents; you may propose a rename, but only as a proposal,
   with the list of call sites it would touch.
4. **Find unused functions.** For every function in `camera_uncertainties.jl` (including those that
   would be moved there) and every function defined in the three scripts, search the whole project for
   uses. For example, `sample_covariance_from_pose_cache` in `uncertainty_camtodrop.jl` appears not to
   be called. Propose for deletion only functions that are not used anywhere. For functions used only
   outside the three scripts, list where they are used.
5. **Find duplicates.** Identify functions that duplicate or nearly duplicate each other or existing
   functions in `camera.jl`. Propose merges, with estimated lines saved and the risk of changing
   results. Propose only; do not change anything.

Report the move plan, the module approach, deletion candidates and merge proposals, then
**wait for my approval.**

## Phase 2 — move the approved functions

- Before any change, copy the current `uncertainty_camtodrop.jld2` to
  `uncertainty_camtodrop_before.jld2` in the same folder.
- Move only the approved functions, without changing their behavior. Update includes, imports, exports
  and calls in the scripts. Delete only functions I explicitly approved for deletion.
- Verify nothing changed:
  - all tests pass;
  - `uncertainty_camtodrop.jl` runs (same seed), and every numeric key in the new
    `uncertainty_camtodrop.jld2` matches `uncertainty_camtodrop_before.jld2` to numerical precision.
    Report the maximum absolute difference per key;
  - `Camera_intrinsics.jl` and `Camera_extrinsics.jl` run and reproduce their saved outputs. If a run
    would take long (e.g. the 500-draw intrinsic bootstrap), ask me first or propose a shorter check.
- Report what was moved, deleted and changed, with the verification results.

## Phase 3 — overview document

Write `docs/uncertainty_camtodrop_overview.md` reflecting the code *after* Phase 2. Reader: me and an AI
assistant I'll share it with, familiar with the physics but not the code. Use LaTeX math, reference
code as `file.jl:line`, and say explicitly where you are unsure. No padding.

1. **Summary** (max half a page): goal, approach, key results with numbers.
2. **Theory and models**, key equations with short explanations, and the assumptions each relies on:
   forward and inverse projection chains (OpenCV distortion [k1, k2, p1, p2, k3], frame composition with
   the relative object pose, meridian-plane intersection); intrinsic calibration and its bootstrap
   covariance; extrinsic pose estimation from ArUco corners, corner fitting and its bootstrap;
   implicit-function-theorem derivative of the undistortion; linearized propagation (camera term and
   pixel term); projection onto contour normals; independent Monte Carlo vs joint bootstrap (what is
   sampled, which noise basis, why the joint version captures intrinsic–extrinsic coupling); eigenmode
   decomposition into rigid-body, scale and shape.
3. **Code structure**: per file, line count and purpose; a table of every function (name, file, purpose
   in one sentence, main inputs/outputs, category: core math / plumbing / test helper); a Mermaid
   data-flow diagram of each of the three scripts showing which functions they use; configurable
   parameters and current values; contents of each saved .jld2 output (keys, shapes, units, frame,
   parameter ordering).
4. **Remaining simplification proposals** that were not carried out, with estimated savings and risk.
5. **Validation**: every test and check, what it verifies, its result; bugs found during development
   (y-residual typo in undistortion, bootstrap noise-basis mismatch), how detected and resolved.
6. **Results**: table of key numbers with the configuration they hold for (zoom5, draws, σ_pix, noise
   basis).
7. **Assumptions, conventions and limitations**, including the rotation perturbation convention (right
   perturbation), virtual contour instead of observed data, σ_pix as placeholder, validity limited to
   the zoom5 calibration, and that `extrinsic_uncertainty_zoom5.jld2` is no longer consistent with this
   analysis.
8. **Status**: what is done, partially done, not done.
9. **Open decisions and next steps**: inspect `main_BI.jl` and report how the relative object pose is
   handled (fixed or inferred, priors) and what that means for the rigid-body modes of the camera
   covariance; outline how the joint covariance could enter the likelihood (full covariance with
   cross-point correlations, low-rank camera part plus pixel noise, Woodbury identity, normal distances
   or (r, z)); other steps (observed contours, measured σ_pix, other zoom settings). List decisions
   needing my input separately from steps you could do yourself.

When finished, reply only with the file path and a three-line summary.