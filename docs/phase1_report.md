# Phase 1 report: uncertainty reorganization

## Scope and status

This report records the analysis and decisions for Phase 1. No source, test, or
data files were changed during the analysis. The only change requested at this
point is creation of this report; Phase 2 remains unapproved.

The scoped files are:

- `ProjectScripts/Camera_intrinsics.jl`
- `ProjectScripts/Camera_extrinsics.jl`
- `ProjectScripts/uncertainty_camtodrop.jl`
- `src/camera_uncertainties.jl`

`src/camera.jl` and all other files are out of scope for code changes.

## Script entry points and test inclusion

| Script | Current entry style and main work |
|---|---|
| `ProjectScripts/Camera_intrinsics.jl` | Top-level code; no named functions and no `main()` guard. Runs subset-calibration analysis, a 500-draw intrinsic bootstrap, a convergence study, and plotting. |
| `ProjectScripts/Camera_extrinsics.jl` | Top-level code; no named functions and no `main()` guard. Detects/fits marker corners, performs the configured extrinsic bootstrap, reports uncertainty, and saves plots/results. |
| `ProjectScripts/uncertainty_camtodrop.jl` | Defines `main()` at line 368 and calls it under the entry guard at line 708. Keep its parse-safe spelling `abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()`. The parentheses around `@__FILE__` are needed to prevent the macro from consuming the `||` expression. |

The user's clarification is that `Camera_intrinsics.jl` and
`Camera_extrinsics.jl` must remain top-level scripts: do not add `main()` or
guards to them. The guard requirement applies only to
`uncertainty_camtodrop.jl`.

A project-wide search found no file that includes any of these three scripts.
In particular, no test includes them. `test/uncertainty_camtodrop.jl` includes
`src/forwardmodel.jl`, not the analysis script, so running tests does not
trigger its full analysis.

## Function classification and proposed moves

`Camera_intrinsics.jl` and `Camera_extrinsics.jl` define no named functions.
Their top-level statements are script orchestration and stay in place.

`uncertainty_camtodrop.jl` defines the following functions:

| Function and location | Classification | Phase 2 disposition |
|---|---|---|
| `require_keys` (`:59`) | Input validation/plumbing | Keep in script. |
| `resample_curve` (`:65`) | Geometric contour preprocessing | Keep in script. |
| `stack_rz` (`:79`), `unstack_rz` (`:89`) | Data layout helpers for contour uncertainty samples | Move to `IntrinsicUncertainty`. |
| `contour_normals` (`:94`) | Computes contour directions used to project covariance | Move. |
| `normal_covariance` (`:110`), `normal_sigmas` (`:122`) | Propagates covariance to normal distances and standard deviations | Move. |
| `block_diagonal_extrinsic` (`:125`) | Assembles a block covariance | Move. |
| `rotation_log` (`:132`) | Converts pose increments to rotation-vector coordinates for covariance analysis | Move. |
| `summarize_pose_samples` (`:142`) | Estimates pose means and a full six-parameter sample covariance | Move. |
| `correlation_matrix` (`:170`) | Converts covariance to correlation form | Move. |
| `sample_covariance_from_pose_cache` (`:176`) | Estimates pose covariance from saved bootstrap samples | Unused; follow the user's deprecation instruction below rather than moving or deleting. |
| `load_or_fit_aruco_corners` (`:210`) | Cache/image I/O and orchestration of detection and corner fits | Keep in script; it depends on script paths/configuration and image-processing packages. |
| `make_normal_band` (`:285`) | Constructs a plotting shape | Keep in script. |
| `apex_radius_measure` (`:294`) | Computes a geometric contour measurement; uncertainty is calculated from its samples elsewhere | Keep in script. |
| `covariance_terms` (`:303`) | Computes linearized camera, pixel, offset, and total covariance terms | Move. |
| `sample_inverse_contour` (`:324`) | Samples camera parameters and inverse-projects pixels into contour coordinates | Move with the required interface change below. |
| `print_leading_modes` (`:336`) | Decomposes covariance into rigid, scale, and shape modes | Move. |
| `main` (`:368`) | End-to-end analysis orchestration, configuration use, output and plotting | Keep in script. |

### Dependencies of moved functions

The moved math helpers use `LinearAlgebra` and/or `Statistics`, already
available to the nested module. `print_leading_modes` also uses `@printf`, so
`Printf` must be imported in `camera_uncertainties.jl`.

Two moved functions refer to script-local names:

- `summarize_pose_samples` uses the script's `UQ` alias to call
  `_right_perturbed_rotation`. In the module, call that module-local function
  directly; no configuration argument is needed.
- `sample_inverse_contour` reads `PLANE_COORDINATE_MM` and uses the `UQ` alias.
  The user approved moving it with a required keyword argument
  `plane_coordinate::Real` (no default), consistent with
  `inverse_project_camtodrop`. Pass `PLANE_COORDINATE_MM` explicitly at every
  script call. Replace the script alias with a module-local call.

No other proposed move reads a script configuration constant. The unused
`sample_covariance_from_pose_cache` reads `EXTRINSIC_DIR` and `ZOOM_LEVEL`, but
will be commented out rather than moved, per the user's instruction.
`load_or_fit_aruco_corners` also uses several script path constants and
`FileIO`/`Images`/`Interpolations`; it stays in the script.

### Ambiguous cases

- `sample_inverse_contour` is clearly uncertainty-related, but its plane
  coordinate is script configuration. The user approved the explicit required
  keyword argument and explicit call-site values.
- `load_or_fit_aruco_corners` obtains corner covariance by calling the
  uncertainty-aware `fit_corner`, but combines that with image detection,
  cache loading/saving, and script-specific paths. Keep this orchestration in
  the script; the reusable fitting implementation already resides in
  `camera_uncertainties.jl`.
- `apex_radius_measure` is applied to Monte Carlo contours, but it computes a
  geometric observable rather than propagating or estimating parameter
  uncertainty. Keep it in the script.
- `resample_curve` prepares the nominal contour; it does not sample uncertain
  parameters. Keep it in the script.
- `make_normal_band` is a plotting helper and stays in the script.

## Module approach

Keep the existing `Camera.IntrinsicUncertainty` module. It is defined in
`src/camera.jl` and includes `src/camera_uncertainties.jl`. The scripts already
access functions by module qualification (`Camera.IntrinsicUncertainty` or the
`UQ` alias), so moved functions need no exports from the nested module.

Do not rename the module. A rename would touch the declaration, export and
include in `src/camera.jl`, plus references in `Camera_intrinsics.jl`,
`Camera_extrinsics.jl`, `uncertainty_camtodrop.jl`, and
`test/uncertainty_camtodrop.jl`. Since `camera.jl` is explicitly out of scope,
retain the current name.

## Project-wide use audit

The audit covered all named functions in the three scripts and all functions
in `camera_uncertainties.jl`, including local named helpers.

### Script-defined functions

All script helpers proposed for moving are used by `main()`. The unused
`sample_covariance_from_pose_cache` has no project call sites. The remaining
script helpers are used in `main()` or by another helper in this script. The
two other scripts have no named functions to audit.

### Functions in `camera_uncertainties.jl`

The following inventory records the project use found for each top-level
function. Uses described as internal are calls within this source file.

| Function(s) | Uses found |
|---|---|
| `_coerce_camera_matrix`, `_coerce_distortion_coeffs` | Internal calibration/bootstrap pipeline and `Camera_intrinsics.jl`. |
| `detect_aruco_corners` | `Camera_extrinsics.jl` and `uncertainty_camtodrop.jl`. |
| `aruco_uncertainty_setup` | No project call sites. |
| `sigmoid` | Internal corner-diagnostic plotting. |
| `fit_corner` | `Camera_extrinsics.jl`, `uncertainty_camtodrop.jl`, and corner-diagnostic helpers. |
| `plot_aruco_corner_corrections`, `plot_aruco_corner_diagnostic`, `plot_aruco_edge_profile_residuals` | `Camera_extrinsics.jl`. |
| `bootstrap_extrinsics`, `extrinsic_covariance` | `Camera_extrinsics.jl`. |
| `list_camera_images`, `random_image_subset` | Internal subset-study pipeline and `Camera_intrinsics.jl`. |
| `calibrate_camera_from_images` | Internal subset-study pipeline and `Camera_intrinsics.jl`. |
| `compute_reprojection_residuals` | Internal reprojection-bootstrap pipeline and `Camera_intrinsics.jl`. |
| `bootstrap_calibration` | Internal reprojection-bootstrap pipeline. |
| `run_subset_study`, `run_reprojection_bootstrap` | `Camera_intrinsics.jl`. |
| `plot_reprojection_residuals` | Internal reprojection-bootstrap pipeline and `Camera_intrinsics.jl`. |
| `summarise_by_subset`, `plot_subset_influence`, `plot_subset_seed_values` | `Camera_intrinsics.jl`. |
| `plot_extrinsic_distributions`, `plot_extrinsic_bootstrap_stability`, `plot_extrinsic_convergence_mean`, `save_extrinsic_uncertainty` | `Camera_extrinsics.jl`. |
| `_camtodrop_primal` (both methods), `_right_perturbed_rotation` | Internal projection/Jacobian helpers; also used by the target script and/or `test/uncertainty_camtodrop.jl`. |
| `inverse_project_camtodrop`, `forward_project_camtodrop`, `camtodrop_jacobians` | `uncertainty_camtodrop.jl`, `test/uncertainty_camtodrop.jl`, and internal Jacobian calls where applicable. |
| `undistortion_implicit_parameter_jacobian` | Only `test/uncertainty_camtodrop.jl`; retain it. |
| `sample_gaussian`, `joint_intrinsic_extrinsic_bootstrap` | `uncertainty_camtodrop.jl`. |

The file also has local named helpers: `make_ellipse` is used within
`plot_aruco_corner_corrections`; `to_pixels` within
`plot_aruco_edge_profile_residuals`; `_fail` within
`calibrate_camera_from_images`; `_panel` within `plot_subset_influence`;
`_seed_panel` within `plot_subset_seed_values`; and `get_result_value` within
`save_extrinsic_uncertainty`. These have no independent project call sites.

### Approved handling of unused functions

The user explicitly instructed that neither unused function be deleted. In
Phase 2:

1. Comment out `sample_covariance_from_pose_cache` in
   `ProjectScripts/uncertainty_camtodrop.jl`.
2. Comment out `aruco_uncertainty_setup` in `src/camera_uncertainties.jl`.
3. Put this exact line directly above each block comment:

   ```julia
   # DEPRECATED: unused as of the uncertainty reorganization; candidate for deletion.
   ```

Use Julia block comments (`#= ... =#`). List both functions in the Phase 2
report and in `docs/uncertainty_camtodrop_overview.md`.

## Duplicate and merge proposals

These are proposals only. No merge is approved or performed in Phase 1.

| Similarity | Proposal | Estimated savings | Risk |
|---|---|---:|---|
| `_coerce_camera_matrix` and `_coerce_distortion_coeffs` are effectively equivalent in `camera.jl` and `camera_uncertainties.jl`. | Reuse the parent `Camera` implementations from the nested module. | About 35–45 lines. | Low; confirm `nothing`, matrix, and nested-vector inputs retain identical behavior. |
| `list_camera_images` and `random_image_subset` are near-identical in both files. | Reuse the exported `Camera` functions and update nested/script call sites. | About 45–55 lines. | Low; verify seeded sampling and extension behavior. |
| `bootstrap_extrinsics` and `joint_intrinsic_extrinsic_bootstrap` repeat parts of the solvePnP setup and draw loop. | Consider extracting only their common low-level solve/geometry code. | About 20–30 lines. | Moderate; preserve each RNG stream, intrinsic sampling, and corner-noise basis exactly. |
| The reprojection-residual and camera-calibration routines in `camera.jl` resemble the nested versions. | No merge proposed. | — | The nested residual routine also returns a 2x2 covariance; nested calibration returns object/image points, poses, image size, and intrinsic standard deviations required by its bootstrap. |
| The two `bootstrap_calibration` routines have similar names and purpose. | Do not merge. | — | High: the nested routine samples the full correlated 2x2 corner covariance; the `camera.jl` routine samples independent x/y noise. A merge could change uncertainty results. |
| `extrinsic_covariance` and `summarize_pose_samples` both estimate pose statistics. | Do not merge. | — | They use different rotation-mean conventions and output contracts; `summarize_pose_samples` retains the full rotation/translation cross-covariance. |

The script's `load_or_fit_aruco_corners` also overlaps in purpose with the
corner-fitting sequence in `Camera_extrinsics.jl`, but that sequence is
top-level rather than a function and uses different configuration/cache
settings. No merge is proposed.

## Phase 2 prerequisites and validation

The existing
`data/test_droplets/camtodrop_uncertainty/uncertainty_camtodrop.jld2` was
confirmed present, so the required pre-change copy can be made before any
Phase 2 edits.

After approval, Phase 2 should move only the approved functions, apply the
required plane-coordinate keyword, preserve the two top-level script entry
styles, and comment out the two unused functions as specified above. Run tests
in the project environment (`--project=.`). Before running the 500-draw
intrinsic or 1000-draw extrinsic scripts, ask if the run is expected to take
long. For the cam-to-drop script, compare every numeric output key against the
backup using the same seed and report the maximum absolute difference per key.
No tests or analyses were run during this analysis-only phase.
