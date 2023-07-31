module InputParser
    using Parsers, DelimitedFiles, Jutul, OrderedCollections, Dates
    export parse_deck_file

    include("deckinput/parser.jl")
end
