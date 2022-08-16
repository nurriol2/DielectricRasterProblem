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

module RectangleAbstraction 

#=
    This module provides methods that are useful for representing pixels and
    Squares as generic rectangles. 

    In addition, an algorithm is implemented to find the intersectional area
    between two generic rectangles. 
=# 


#=
    Compute the lines forming the boundary of a (Square) rectangle from a Square

    Ex: Top surface -> y = (Square center y coordinate) + (1/2 Square width)
=#
right_edge(s::Square) = s.x + s.w/2
left_edge(s::Square) = s.x - s.w/2
top_edge(s::Square) = s.y + s.w/2
bot_edge(s::Square) = s.y - s.w/2

# Rectangle representation for the pixels found in `GriddedArray.pixels`
function right_edge(A::GriddedArray, index::Int64)
    # Find x coord pixel center
    xcenter, _ = pixel_center(A, index)
    # Find pixel width
    w = pixel_width(A)
    # Compute edge loc
    xcenter + w/2
end

function left_edge(A::GriddedArray, index::Int64)
    # Find x coord pixel center
    xcenter, _ = pixel_center(A, index)
    # Find pixel width
    w = pixel_width(A)
    # Compute edge loc
    xcenter - w/2
end

function top_edge(A::GriddedArray, index::Int64)
    # Find y coord pixel center
    _, ycenter = pixel_center(A, index)
    # Find pixel height
    h = pixel_height(A)
    # Compute edge loc
    ycenter + h/2
end

function bot_edge(A::GriddedArray, index::Int64)
    # Find y coord pixel center
    _, ycenter = pixel_center(A, index)
    # Find pixel height
    h = pixel_height(A)
    # Compute edge loc
    ycenter - h/2
end

# Algorithm for calculating intersectional area between generic rectangles
function w_prime(A::GriddedArray, index::Int64, s::Square)
   min(right_edge(s), right_edge(A, index)) - max(left_edge(s), left_edge(A, index)) 
end

function h_prime(A::GriddedArray, index::Int64, s::Square)
   min(top_edge(s), top_edge(A, index)) - max(bot_edge(s), bot_edge(A, index))
end

function overlapping_area(A::GriddedArray, index::Int64, s::Square)
    max(0, w_prime(A, index, s)) * max(0, h_prime(A, index, s))
end

#= 
    For future problems, it will be useful to construct the raster matrix without
    without modifying the underlying matrix.
=#
function raster_area(A::GriddedArray, s::Square)
    
    # Create a copy of the unmodified GriddedArray
    # FIXME:  Potentially expensive operation
    A′ = deepcopy(A)
    
    # Calculate the total area of a single pixel
    pixel_area = total_pixel_area(A′)

    # Iterate over the pixels in A
    for index in eachindex(A′.pixels)
        # Calculate the overlapping area
        overlap_area = overlapping_area(A′, index, s)
        # Calculate the target
        # Assign the value of the target to the current element in A
        A′.pixels[index] = (overlap_area / pixel_area)
    end
    A′
end


end # RectangleAbstraction

## Q1
function raster_area!(A::GriddedArray, s::Square)
    # Calculate the total area of a single pixel
    pixel_area = total_pixel_area(A)

    # Iterate over the pixels in A
    for index in eachindex(A.pixels)
        # Calculate the overlapping area
        overlap_area = overlapping_area(A, index, s)
        # Assign the ratio of areas to the current pixel in A
        A.pixels[index] = (overlap_area / pixel_area)
    end
    A
end


## Q2
function raster_scalar_dielectric!(A::GriddedArray, s::Square)
    
    #TODO: Mutate the input matrix `A`  

    # Intersectional area matrix
    O = raster_area(A, s).pixels 
    # "Compliment" to the intersectional area matrix
    C = ones(size(A.pixels)) - O 
    # Introduces the dielectric constant in matrix form for formula below 
    K_square = s.dielectric .* ones(size(A.pixels))

    # This formula computes the solution to the problem
    # (A.pixels .* C) + (K_square .* O)
    
    # An extra variable is created for the operation below
    solution_matrix = (A.pixels .* C) + (K_square .* O)

    #=
        HACK:  The only way I've found to actually mutate a matrix in place 
                is by iterating over its elements and setting the value.
                This approach can be extremely slow.
                With more practice writing Julia, I'm confident there is an
                idiomatic way to do this common operation. 
                The built-in `replace!()` method looked promising.
    =#
    for index in eachindex(solution_matrix)
        A[index] = solution_matrix[index]
    end
    A
end

## Q3
function raster_harmonic!(A::GriddedArray, s::Square)
    # Intersectional area matrix
    O = raster_area(A, s).pixels 
    # "Compliment" to the intersectional area matrix
    C = ones(size(A.pixels)) - O 
    # Introduces the dielectric constant in matrix form for formula below 
    K_square = s.dielectric .* ones(size(A.pixels))

    # This solution does not modify the matrix
    ((A.pixels) .* K_square) ./ ((K_square .* C) .+ (A.pixels .* O))
    
    #=
        HACK:  The only way I've found to actually mutate a matrix in place 
                is by iterating over its elements and setting the value.
                This approach can be extremely slow.
                With more practice writing Julia, I'm confident there is an
                idiomatic way to do this common operation. 
                The built-in `replace!()` method looked promising.
    =#
    
end

module FuzzyComparisons

function fuzzy_eq(a::Tuple, b::Tuple)
    (a[1] ≈ b[1]) && (a[2] ≈ (b[2]))
end

function fuzzy_leq(left, right)
    (left < right) || (left ≈ right)
end

function fuzzy_in_open_range(lep, x, rep)

    #=
        Check if x is inside the open range (lep, rep).
        
        Examples:
            ```
            julia> fuzzy_in_open_range(1.0, 1.5, 2.0)
            true

            julia> fuzzy_in_open_range(1.0, 1.0, 2.0)
            false

            julia> fuzzy_in_open_range(1.0, 2.0, 2.0)
            false

            julia> fuzzy_in_open_range(2.0, 1.5, 1.0)
            ERROR:  ArgumentError
        ```

        Raises:
            ArgumentError: If endpoints are equal OR `rep` < `lep`
    =#
    
    # Since the range is open, x cannot be lep or rep
    
    # TODO:  Write custom errors with explanation of argument order
    if fuzzy_leq(rep, lep)
        throw(ArgumentError)
    end
   
    # Check if the point is one of the endpoints
    if x ≈ lep
        return false
    end
    
    if x ≈ rep
        return false
    end
       
    # Only eval to true when x is strictly in the open range
    fuzzy_leq(lep, x) && fuzzy_leq(x, rep)
    
end

end # FuzzyComparisons

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
