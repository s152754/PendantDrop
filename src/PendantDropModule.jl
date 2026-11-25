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
    using Dierckx
    using Contour
    using Interpolations
    using LsqFit
    using Statistics
    using LinearAlgebra

    include("data_processing.jl")
    # include all other function files as needed for the notebook

    # exports from "data_processing.jl"
    export test, load_file, RectangleSelector, extract_roi, find_max_gradients, build_binary_mask, extract_boundary_coordinates, trim_needle, fit_spline

end # module PendantDropModule
