using Test
import JutulDarcy.InputParser: clean_include_path, parse_defaulted_line
@testset "InputParser" begin
    @test clean_include_path("", " MYFILE") == "MYFILE"
    @test clean_include_path("/some/path", " MYFILE") == joinpath("/some/path", "MYFILE")
    @test clean_include_path("/some/path", " 'MYFILE'") == joinpath("/some/path", "MYFILE")
    @test clean_include_path("/some/path", " ./MYFILE") == joinpath("/some/path", "MYFILE")
    @test clean_include_path("/some/path", " './MYFILE'") == joinpath("/some/path", "MYFILE")
    @test clean_include_path("/some/path", " 'INCLUDE/file.txt' /  (Some comment)") == "/some/path/INCLUDE/file.txt"

    @test parse_defaulted_line("3.0 2* 7", [1.0, 2, 3, 4]) == [3.0, 2, 3, 7]
    @test parse_defaulted_line("2.0", [1.0, 2, 3, 4]) == [2.0, 2, 3, 4]
    @test parse_defaulted_line("5 *", [1, 2]) == [5, 2]
    @test parse_defaulted_line("5  2*", [1, 2, 3]) == [5, 2, 3]
    @test parse_defaulted_line("5, 2", [1, 2]) == [5, 2]
    @test parse_defaulted_line("5   HEI", [1, "Hallo"]) == [5, "HEI"]
    @test parse_defaulted_line("2*", [1, "Hallo"]) == [1, "Hallo"]
end
