# Task: hand-over point for the student's edge detection + total data uncertainty

A student has her own script that detects the drop edge and estimates its uncertainty. She will deliver
her output as a JLD2 FILE; do not integrate her code. Prepare the pipeline to read that file and
combine her edge uncertainty with the camera uncertainty into a total covariance, and write a README for
her. Credits are limited: work efficiently, do not run long analyses, do not compile LaTeX, and do not
stop for approval; where something is ambiguous, choose the simplest option and document it in the
README.

## 1. File contract (her output): one `.jld2` file with these keys
Required:
- `"u"`, `"v"`: Vector{Float64}, pixel coordinates of the edge points. Determine the exact convention
  the pipeline uses (1-based Julia pixels, position of the centre of pixel (1,1), raw distorted image)
  from the code, e.g. the +1 shift applied to OpenCV corners and the output of
  `forward_project_camtodrop`; cite file:line.
- `"sigma_n_px"`: Vector{Float64}, standard deviation of the edge position normal to the edge, in
  pixels, one value per point.
- `"image_file"`: String, the image the contour was detected in.
Optional correlation between points (at most one of these):
- `"corr_length_px"`: Float64, correlation length along the edge for an exponential model
  ρ_ij = exp(−|s_i − s_j| / ℓ), with s the arc length along the detected edge in pixels;
- `"corr_matrix"`: Matrix{Float64}, the full N×N correlation matrix of the normal edge errors.
If neither is given, edge errors are uncorrelated.
Also specify the expected point ordering and which part of the contour (one side or both sides; apex to
needle or reverse), matching what the downstream code expects from `virtual_u`/`virtual_v`. The
pipeline must ignore extra keys, so she can store her own metadata in the same file.

## 2. Reading and combining (computation in `src/camera_uncertainties.jl`, in the fitting sections;
##    configuration in `ProjectScripts/uncertainty_camtodrop.jl`)
- Add `CONTOUR_SOURCE = :virtual` (alternative `:observed`) and `EDGE_DATA_FILE` to the script.
- `:virtual` must give exactly the current results (back up the .jld2 and compare all keys).
- `:observed`: read and validate the file; build the edge covariance in pixels per point as
  σ_n² · n_img n_imgᵀ (normal-direction only, with n_img the image-space edge normal computed from the
  points), with off-diagonal blocks from the correlation if given; propagate it to the model plane with
  the existing pixel Jacobian J_p; compute the joint camera covariance at her pixel locations; project
  both onto the model-plane contour normals; and form the total normal covariance
  Σ_total = Σ_camera + Σ_edge (N×N, mm²), keeping all correlations.
- Save the contour in mm (r, z), Σ_camera, Σ_edge and Σ_total (normal direction, mm²), plus per-point σ
  for each, under new keys.
- Where a downstream step assumes the virtual contour (comparisons with `virtual_r`/`virtual_z`, the
  round-trip check), skip it in `:observed` mode with a printed note.
- Add a figure: σ per point along the edge for camera, edge and total; and the correlation matrix of
  Σ_total.

## 3. Tests
- Contract test for an edge-data file: required keys and types, equal lengths, finite values, points
  inside the image, expected ordering, σ positive, correlation (if given) valid: `corr_length_px`
  positive, or `corr_matrix` N×N, symmetric, unit diagonal, positive semidefinite.
- Plumbing test: an edge file made from the virtual pixels with a constant `sigma_n_px` must reproduce
  the existing camera covariance terms, and an edge term consistent with the current isotropic pixel
  term projected onto the normals.
- An example file `docs/example_edge_data.jld2` generated from the virtual contour, valid under the
  contract, plus the few lines of Julia that created it (in the README), so she can copy the pattern.
- Confirm all existing tests pass.

## 4. README for the student: `docs/README_edge_uncertainty.md` (max ~2 pages)
- What to deliver: the file contract, with the pixel convention and ordering explained plainly and why
  they matter (1 px ≈ 8 µm on the drop; a convention mismatch shifts everything without an error).
- A short Julia example of writing a valid file with `jldsave`.
- How to check her file with the contract test, and how to run the pipeline in `:observed` mode.
- What happens to her σ: converted to mm via the camera model (camera parameters plus the assumption
  that the drop lies in the marker plane), combined with the camera covariance, and what the resulting
  total covariance and its correlations mean (camera: correlated along the whole edge; edge: locally
  correlated).
- How to provide correlation along the edge and why it matters (ignoring it makes the final
  uncertainty too optimistic).
- Where the results and figures appear, and which documents explain the rest.

When finished, report briefly: the pixel convention and ordering you determined (with file:line), the
verification results, and the list of files changed (for committing).