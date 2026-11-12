# Keep this list in sync with the test matrix in .github/workflows/CI.yml.
const CI_TEST_GROUPS = (
    "core", "properties", "flow", "wells", "well_controls", "backends",
    "solvers", "fractures", "partitioned", "thermal", "thermal_input",
    "compositional", "adjoint_parameters", "adjoint_controls",
)
const CI_TEST_TAGS = Symbol.("ci_" .* CI_TEST_GROUPS)
const CI_TEST_GROUP = get(ENV, "JUTULDARCY_TEST_GROUP", "all")
const MANUAL_TEST_TAGS = Symbol.(filter(!isempty, strip.(split(get(ENV, "JUTULDARCY_TEST_TAGS", ""), ','))))

if CI_TEST_GROUP != "all" && !(CI_TEST_GROUP in CI_TEST_GROUPS)
    error("Unknown JUTULDARCY_TEST_GROUP=$(repr(CI_TEST_GROUP)). Choose all or one of $(CI_TEST_GROUPS).")
end

function testitem_in_ci_group(item)
    groups = filter(tag -> startswith(string(tag), "ci_"), item.tags)
    if length(groups) != 1 || !(only(groups) in CI_TEST_TAGS)
        error("Test item $(repr(item.name)) must have exactly one valid CI group tag, got $(groups).")
    end
    in_group = CI_TEST_GROUP == "all" || only(groups) == Symbol("ci_", CI_TEST_GROUP)
    return in_group && all(tag -> tag in item.tags, MANUAL_TEST_TAGS)
end
