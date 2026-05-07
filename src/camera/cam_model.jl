## determine intrinsic matrix K
function getK_fn(dir_cam, file_cam)

    file_out = joinpath(dir_cam, file_cam)
    @load cam_matrix

    return cam_matrix
end

## determine rotation matrix R| t

## translate model coordinates to image coordinates