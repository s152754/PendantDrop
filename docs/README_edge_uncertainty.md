# Edge uncertainty hand-over

The pipeline reads your measured contour and per-point edge uncertainty from a
JLD2 file; it does not run or integrate your edge detector. The `:observed`
mode combines your edge uncertainty with the joint camera-calibration and
ArUco-pose uncertainty.

## What to deliver

Deliver one `.jld2` file with:

| Key | Required value |
|---|---|
| `u`, `v` | `Vector{Float64}` coordinates in pixels |
| `sigma_n_px` | `Vector{Float64}` positive 1σ standard deviations normal to the edge, in px |
| `sigma_description::String` | Non-empty description of included error sources and how `sigma_n_px` was estimated |
| `image_file` | `String` path to the full-frame drop image |

Using one side is the current scope; using both sides is a possible later
extension.

The vectors must have the same length and contain finite values. Points must
be inside the image. `sigma_n_px` is a 1σ standard deviation; which error
sources it includes is your choice, but report those sources and how it was
estimated in `sigma_description`. The image must exist and have the same
dimensions as the ArUco image. Use a full-frame drop image taken with the
zoom5 setup (same camera, zoom, and focus), with neither the camera nor marker
moved since the ArUco image. Any file name is fine; relative paths are resolved
from the repository root. Do not crop, resize, or undistort the image before
measuring the contour.

**Pixel convention matters.** Store Julia-style 1-based image coordinates:
`u` is the column, `v` is the row, and the centre of array pixel `(1,1)` is
coordinate `(1,1)`. A one-pixel convention error shifts the whole contour
without necessarily raising an error; 1 px is about 8 µm on the drop. The
ArUco loader adds 1 to OpenCV corners before sampling Julia image arrays
(`src/camera_uncertainties.jl:709-711`). In contrast,
`forward_project_camtodrop` returns OpenCV-style zero-based coordinates from
the camera intrinsics (`src/camera_uncertainties.jl:1116,1140`;
`src/camera.jl:966-967`). The observed pipeline subtracts 1 from your stored
coordinates before using that projection model (`ProjectScripts/uncertainty_camtodrop.jl:221`).
Coordinates must refer to pixel centres. If your detector returns 0-based
coordinates (e.g. Python/OpenCV), add 1 to both u and v.

Provide the **RIGHT-hand side of the drop silhouette as seen in the image**,
ordered **apex to needle**, and exclude the needle edges. The virtual
contour follows the `PD_fn` integration order: it starts near `r=0`, progresses
along the drop, and ends at the needle-radius crossing
(`src/forwardmodel.jl:135,247-249`; `ProjectScripts/uncertainty_camtodrop.jl:192`).
The radius can increase and then decrease, so do not sort points by `r`.
The pipeline accepts either sign of the projected meridian radius and checks
the endpoints and one-sided ordering.

Optionally provide exactly one correlation key:

* `corr_length_px::Float64` > 0: an exponential model,
  `rho[i,j] = exp(-abs(s[i]-s[j]) / corr_length_px)`, where `s` is cumulative
  arc length along your ordered image contour.
* `corr_matrix::Matrix{Float64}`: the full correlation matrix. It must be
  finite, symmetric, unit-diagonal, and positive semidefinite.

If neither key is present, the pipeline assumes uncorrelated edge errors.
Extra keys are ignored, so you may include your own metadata.

## Example file

The checked-in [example_edge_data.jld2](./example_edge_data.jld2) uses the
virtual contour as a stand-in and a constant 0.3 px uncertainty. These few
lines generated it:

```julia
using JLD2
reference = load("data/test_droplets/camtodrop_uncertainty/uncertainty_camtodrop.jld2")
u, v = reference["virtual_u"] .+ 1, reference["virtual_v"] .+ 1
jldsave("docs/example_edge_data.jld2"; u, v,
        sigma_n_px=fill(0.3, length(u)),
        sigma_description="Constant 0.3 px 1σ, assumed for this example.",
        image_file="data/test_droplets/zoom5.JPG", corr_length_px=8.0)
```

The `+1` converts the virtual OpenCV coordinates to the file contract's
1-based Julia coordinates.

## Check and run

Run the focused contract and covariance tests:

```powershell
julia --project=. .\test\uncertainty_camtodrop.jl
```

To validate your file, including its projected apex-to-needle ordering, run
from the repository root:

```powershell
julia --project=. .\test\check_edge_uncertainty_file.jl .\path\to\your_edge_data.jld2
```

The checker uses the saved virtual-run pose in
`data/test_droplets/camtodrop_uncertainty/uncertainty_camtodrop.jld2`.
Then set `CONTOUR_SOURCE = :observed` and `EDGE_DATA_FILE` near the top of
`ProjectScripts/uncertainty_camtodrop.jl`, and run:

```powershell
julia --project=. .\ProjectScripts\uncertainty_camtodrop.jl
```

Your `sigma_n_px` is converted to model-plane uncertainty through the local
pixel Jacobian, camera intrinsics, and the assumption that the drop lies in
the marker plane (`plane_coordinate = 0`). The camera term is estimated at
your pixel locations using the joint intrinsic/pose bootstrap. The total
normal covariance is `Sigma_total = Sigma_camera + Sigma_edge`; its diagonal
gives per-point variance and its off-diagonal entries describe correlations
between contour points.

Camera errors are shared across the contour, so camera-induced correlations
can extend along its full length. Edge errors follow the correlation model
you provide and are often more local. Ignoring real positive edge correlation
can make the final uncertainty too optimistic.

The observed result is written to
`data/test_droplets/camtodrop_uncertainty/uncertainty_camtodrop_observed.jld2`.
The sigma-components figure and total-correlation matrix are saved in the
same directory as `observed_edge_sigma_components.png` and
`observed_edge_total_correlation.png`. For the camera model and covariance
details, see [the cam-to-drop overview](./uncertainty_camtodrop_overview.md)
and [the visual overview](./camera_uncertainty_visual_overview.tex).
