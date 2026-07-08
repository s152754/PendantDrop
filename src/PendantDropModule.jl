"""
PendantDropModule

A Julia package for pendant drop tensiometry data processing and analysis.
"""

module PendantDropModule

    using Base64
    using FileIO
    using HypertextLiteral
    using Images
    using ImageShow
    using JLD2
    using Plots
    using PlutoUI
    using Printf
    using Revise

    include("data_processing.jl")
    # include all other function files as needed for the notebook

    include("testfunctions.jl")

    # exports from "data_processing.jl"
    export test, load_file, RectangleSelector, extract_roi, find_max_gradients, build_binary_mask, extract_boundary_coordinates

    export testing_fn

end # module PendantDropModule
