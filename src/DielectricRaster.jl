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


module SquareOperations
#=
    This module provides methods for querying the relationship between `Square`s
    and points in the plane. 

    For example, determining if a point is on the boundary of a given `Square`
    or if a point is also a corner.


    Throughout, a *surface* is an open range along the boundary of a `Square`.
    
    For example, a `Square` centered at (0, 0) with width 2 (units) has a top 
    surface comprised of all points with y-coordinate = 1 and x-coordinate 
    satisfying -1 < x < 1.
=#

# Canonical distance formula is used througout this module
distance_formula(a::Tuple, b::Tuple) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)

# Corner shorthand:  [(U)pper/(L)ower] [(L)eft/(R)ight] (C)orner
urc(s::Square) = (s.x + s.w/2, s.y + s.w/2)
ulc(s::Square) = (s.x - s.w/2, s.y + s.w/2)
llc(s::Square) = (s.x - s.w/2, s.y - s.w/2)
lrc(s::Square) = (s.x + s.w/2, s.y - s.w/2)
corners(s::Square) = (urc(s), ulc(s), llc(s), lrc(s))


function point_is_corner(p::Tuple, s::Square)
    # Assume point is not a corner
    result = false
    
    # Compare point with coordinates of the corners
    for corner in (urc(s), ulc(s), llc(s), lrc(s))
        if fuzzy_eq(p, corner)
            result = true
        end
    end

    result
end


function which_corner(p::Tuple, s::Square)
    fwd = Dict(
        "urc" => (s.x + s.w/2, s.x + s.w/2),
        "ulc" => (s.x - s.w/2, s.x + s.w/2),
        "llc" => (s.x - s.w/2, s.x - s.w/2),
        "lrc" => (s.x + s.w/2, s.x - s.w/2)
    )
    
    bkd = Dict(v => k for (k, v) in fwd)
    
    bkd[p]
end

# Check if a point is or is not on a specific surface
point_on_top(p::Tuple, s::Square) = (p[2] ≈ (s.y + s.w/2)) && fuzzy_in_open_range(ulc(s)[1], p[1], urc(s)[1])
point_on_bot(p::Tuple, s::Square) = (p[2] ≈ (s.y - s.w/2)) && fuzzy_in_open_range(llc(s)[1], p[1], lrc(s)[1])
point_on_left(p::Tuple, s::Square) = (p[1] ≈ (s.x - s.w/2)) && fuzzy_in_open_range(llc(s)[2], p[2], ulc(s)[2])
point_on_right(p::Tuple, s::Square) = (p[1] ≈ (s.x + s.w/2)) && fuzzy_in_open_range(lrc(s)[2], p[2], urc(s)[2])

# Check if a point is on any surface
function point_is_on_surface(p::Tuple, s::Square)
    result = false
    
    if point_on_top(p, s)
        result = true
    end
        
    if point_on_bot(p, s)
        result = true
    end
    
    if point_on_left(p, s)
        result = true
    end
    
    if point_on_right(p, s)
        result = true
    end
    
    result 
end

function which_surface(p::Tuple, s::Square)
    surface = ""
    
    if point_on_top(p, s)
        surface = "top"
    end
        
    if point_on_bot(p, s)
        surface = "bot"
    end
    
    if point_on_left(p, s)
        surface = "left"
    end
    
    if point_on_right(p, s)
        surface = "right"
    end
    
    surface 
end

# Check if a point is within the boundary of `Square`
function square_contains_point(p::Tuple, s::Square)
    fuzzy_in_open_range(ulc(s)[1], p[1], urc(s)[1]) && fuzzy_in_open_range(llc(s)[2], p[2], ulc(s)[2])
end

#=
    Classifying a point as *inside* a `Square` yields some advantages for
    calculating the closest surface normal.

    Throughout, a point is defined as inside when one of the following statements
    is true

    1. The point is itself the center of `Square`
    2. The point is one of the corners of `Square`
    3. The point lies on one of the `Square` surfaces
    4. The point lies on the interior of the `Square` boundary

=#
function point_is_inside(p::Tuple, s::Square)
    # Assume point is not inside
    result = false
    
    # Check if p is the Square center
    if fuzzy_eq(p, (s.x, s.y))
        result = true
    end
    
    # Check if p is a corner
    if point_is_corner(p, s)
        result = true
    end
    
    # Check if p lies on a surface
    if point_is_on_surface(p, s)
        result = true
    end
    
    # Check if s contains p
    if square_contains_point(p, s)
        result = true
    end
    
    result 
end

# Check if dropping a perpendicular from p intersects with a surface
drop_perp_to_top(p::Tuple, s::Square) = point_on_top((p[1], (s.y + s.w/2)), s)
drop_perp_to_bot(p::Tuple, s::Square) = point_on_bot((p[1], (s.y - s.w/2)), s)
drop_perp_to_left(p::Tuple, s::Square) = point_on_left(((s.x - s.w/2), p[2]), s)
drop_perp_to_right(p::Tuple, s::Square) = point_on_right(((s.x + s.w/2), p[2]), s)

function distance_to_surface(p::Tuple, s::Square, surface)
    surfaces = Dict(
        "top" => (p[1], (s.y + s.w/2)),
        "bot" => (p[1], (s.y - s.w/2)),
        "left" => ((s.x - s.w/2), p[2]),
        "right" => ((s.x + s.w/2), p[2])
    )
    
    if surface == "top"
        # A surface is only a candidate solution if dropping a perpendicular from p lands on the surface
        if drop_perp_to_top(p, s)
            return distance_formula(p, surfaces[surface])
        # This entry will be eliminated when distances are filtered 
        else
            return Inf
        end
    end
    
    if surface == "bot"
        # A surface is only a candidate solution if dropping a perpendicular from p lands on the surface
        if drop_perp_to_bot(p, s)
            return distance_formula(p, surfaces[surface])
        # This entry will be eliminated when distances are filtered 
        else
            return Inf
        end
    end
        
    if surface == "left"
        # A surface is only a candidate solution if dropping a perpendicular from p lands on the surface
        if drop_perp_to_left(p, s)
            return distance_formula(p, surfaces[surface])
        # This entry will be eliminated when distances are filtered 
        else
            return Inf
        end
    end
                            
    if surface == "right"
        # A surface is only a candidate solution if dropping a perpendicular from p lands on the surface
        if drop_perp_to_right(p, s)
            return distance_formula(p, surfaces[surface])
        # This entry will be eliminated when distances are filtered 
        else
            return Inf  
        end
    end

end

function distance_to_corner(p::Tuple, s::Square, corner)
    corners = Dict(
        "ulc" => ulc(s),
        "urc" => urc(s),
        "llc" => llc(s),
        "lrc" => lrc(s)
    )
    distance_formula(p, corners[corner])
end

end # SquareOperations

## Q4
function find_nearest_surface_normal(s::Square, coord)

    normalize!(vector) = vector ./ distance_formula((0, 0), (vector[1], vector[2]))

    solution_map = Dict(
    "top" => [0.0 1.0],
    "bot" => [0.0 -1.0],
    "left" => [-1.0 0.0],
    "right" => [1.0 0.0],
    "urc" => [1/sqrt(2) 1/sqrt(2)],
    "ulc" => [-1/sqrt(2) 1/sqrt(2)],
    "llc" => [-1/sqrt(2) -1/sqrt(2)],
    "lrc" => [1/sqrt(2) -1/sqrt(2)],
    "zero" => [0.0 0.0]
)
    
    # When coord is inside, there are some special cases that immediately give the solution
    if point_is_inside(coord, s)
        # Check if coord is the Square center
        if fuzzy_eq(coord, (s.x, s.y))
            return solution_map["zero"]
        end

        # Check if coord is a corner
        if point_is_corner(coord, s)
            return solution_map[which_corner(coord, s)]
        end

        # Check if coord lies on a surface
        if point_is_on_surface(coord, s)
            return solution_map[which_surface(coord, s)]
        end
    end
    
    # None of the special cases occurred so, iterate through the distances to each target
    
    # Container for candidate targets
    candidate = Dict()
    min_dist = Inf
    min_target = ""
    
    for surface in ("top", "bot", "left", "right")
        dist2surf = distance_to_surface(coord, s, surface)
        # If dist_to_surface ≈ min_dist
        if dist2surf ≈ min_dist
            # Record a tie in candidate
            candidate[surface] = dist2surf
        # If dist_to_surface < min_dist
        elseif dist2surf < min_dist
            # There is a new minimum distance and corresponding target
            for key in keys(candidate)
                delete!(candidate, key)
            end
            # Record the new min_dist
            min_dist = dist2surf
            # Record the new min_target
            min_target = surface
            # Record the best candidate
            candidate[min_target] = min_dist
        end
    end
    

    for corner in ("ulc", "urc", "llc", "lrc")
        dist2cor = distance_to_corner(coord, s, corner)
        # If dist_to_corner ≈ min_dist
        if dist2cor ≈ min_dist
            # Record a tie
            candidate[min_target] = min_dist
        # If dist_to_corner < min_dist
        elseif dist2cor < min_dist
            # There is a new minimum distance and corresponding target
            for key in keys(candidate)
                delete!(candidate, key)
            end
            # Record the new min_dist
            min_dist = dist2cor
            # Record the new min_target
            min_target = corner
            # Record the best candidate
            candidate[min_target] = min_dist
        end
    end
    
    # Average the best candidates
    solution = [0 0]
    for key in keys(candidate)
        solution = (solution .+ solution_map[key])
    end
    solution = solution ./ size(collect(keys(candidate)), 1)
    normalize!(solution)
    
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
