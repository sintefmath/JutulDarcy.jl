using Test
import JutulDarcy: determine_cell_overlap_inside_line
@testset "corner point pillar point overlap" begin
    for start in -10:25
        for increment in 1:25
            top = start
            mid = start + increment
            bottom = mid + increment
            far = bottom + increment
            # A_CONTAINS_B and B_CONTAINS_A
            # A B
            # | |   |
            # | | = |
            # | |   |
            @test determine_cell_overlap_inside_line(top, bottom, top, bottom) == (JutulDarcy.AB_RANGES_MATCH, top:bottom, top:bottom, top:bottom)
            # AB_OVERLAP_A_FIRST and AB_OVERLAP_B_FIRST
            # A B
            #   |   x
            # | | = |
            # | |   x
            @test determine_cell_overlap_inside_line(mid, bottom, top, mid) == (JutulDarcy.AB_OVERLAP_B_FIRST, mid:mid, mid:bottom, top:mid)
            # Reversed case.
            @test determine_cell_overlap_inside_line(top, mid, mid, bottom) == (JutulDarcy.AB_OVERLAP_A_FIRST, mid:mid, top:mid, mid:bottom)
            # TOP_MATCHES_A_LONG and TOP_MATCHES_B_LONG
            # A B
            # | |   |
            # | | = |
            # |     x
            @test determine_cell_overlap_inside_line(top, bottom, top, mid) == (JutulDarcy.TOP_MATCHES_A_LONG, top:mid, top:bottom, top:mid)
            # Reversed case.
            @test determine_cell_overlap_inside_line(top, mid, top, bottom) == (JutulDarcy.TOP_MATCHES_B_LONG, top:mid, top:mid, top:bottom)
            # BOTTOM_MATCHES_A_LONG and BOTTOM_MATCHES_B_LONG
            # A B
            # | x   x
            # | | = |
            # | |   |
            @test determine_cell_overlap_inside_line(top, bottom, mid, bottom) == (JutulDarcy.BOTTOM_MATCHES_A_LONG, mid:bottom, top:bottom, mid:bottom)
            # Reversed case.
            @test determine_cell_overlap_inside_line(mid, bottom, top, bottom) == (JutulDarcy.BOTTOM_MATCHES_B_LONG, mid:bottom, mid:bottom, top:bottom)
            # DISTINCT_A_ABOVE and DISTINCT_B_ABOVE
            # A B
            # |     x
            # |     x
            #   | = x
            #   |   x
            @test determine_cell_overlap_inside_line(top, mid, bottom, far) == (JutulDarcy.DISTINCT_A_ABOVE, 0:-1, top:mid, bottom:far)
            # Reversed case.
            @test determine_cell_overlap_inside_line(bottom, far, top, mid) == (JutulDarcy.DISTINCT_B_ABOVE, 0:-1, bottom:far, top:mid)
        end
    end
end
