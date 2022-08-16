module DielectricRaster

using StaticArrays

struct Square{T}
    x::Float64
    y::Float64
    w::Float64
    dielectric::T
end

struct GriddedArray{T}
    xlims::Tuple{Float64, Float64}
    ylims::Tuple{Float64, Float64}
    pixels::Matrix{T}
end

const TensorGriddedArray = GriddedArray{SMatrix{2, 2, Float64}}

module PixelOperations
#=
    This module provides getters for  pixel characteristics, such as height
    and area.

    These pixels are defined in `GriddedArray.pixels`
=#

function pixel_width(A::GriddedArray)
    ncols = size(A.pixels, 2)
    (A.xlims[2] - A.xlims[1]) / (ncols - 1)

end

function pixel_height(A::GriddedArray)
    nrows = size(A.pixels, 1)
    (A.ylims[2] - A.ylims[1]) / (nrows-1)
end

function pixel_center(A::GriddedArray, index::Int64) 

    nrows = size(A.pixels, 1)
    ncols = size(A.pixels, 2)
    
    Δx = pixel_width(A)
    Δy = pixel_height(A)
    
    xcenter = A.xlims[1] + Δx * ((index - 1.0) % ncols) 
    ycenter = A.ylims[1] + Δy * ((index - 1.0) ÷ ncols)
    
    (xcenter, ycenter)
end

total_pixel_area(A::GriddedArray) = pixel_width(A) * pixel_height(A)

end # PixelOperations

## Q1
function raster_area!(A::GriddedArray, s::Square)
    # implement me
end

## Q2
function raster_scalar_dielectric!(A::GriddedArray, s::Square)
    # implement me
end

## Q3
function raster_harmonic!(A::GriddedArray, s::Square)
    # implement me
end

## Q4
function find_nearest_surface_normal(s::Square, coord)
    # implement me
end

## Q5
function orient_along_normal(A, normal)
    # implement me
    # return `A` oriented in a new reference frame with one axes parallel to `normal`
end

## Q6, OPTIONAL
function raster_tensor_dielectric!(A::GriddedArray{<:SMatrix}, B::GriddedArray, s::Square)
    # implement me
end

export Square, GriddedArray, TensorGriddedArray
export raster_area!, raster_scalar_dielectric!, raster_harmonic!
export find_nearest_surface_normal, orient_along_normal, raster_tensor_dielectric!

end
