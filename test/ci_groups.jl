# Keep these tags in sync with the test matrix in .github/workflows/CI.yml.
const CI_TEST_GROUPS = (
    :core, :flow, :wells, :backends, :solvers, :partitioned, :thermal,
    :compositional, :adjoint_parameters, :adjoint_controls,
)
const CI_TEST_GROUP = get(ENV, "JUTULDARCY_TEST_GROUP", "all")
const CI_TEST_TAGS = Symbol.("ci_" .* string.(CI_TEST_GROUPS))
const MANUAL_TEST_TAGS = Symbol.(filter(!isempty, strip.(split(get(ENV, "JUTULDARCY_TEST_TAGS", ""), ','))))
if CI_TEST_GROUP != "all" && !(Symbol(CI_TEST_GROUP) in CI_TEST_GROUPS)
    error("Unknown JUTULDARCY_TEST_GROUP=$(repr(CI_TEST_GROUP)). Choose all or one of $(CI_TEST_GROUPS).")
end

function testitem_in_ci_group(ti)
    # Check excluded items too, so a missing tag cannot silently drop tests from CI.
    groups = filter(tag -> tag in CI_TEST_TAGS, ti.tags)
    length(groups) == 1 || error("Test item $(repr(ti.name)) must have exactly one CI group tag, got $(groups).")
    in_group = CI_TEST_GROUP == "all" || Symbol("ci_", CI_TEST_GROUP) in groups
    return in_group && all(tag -> tag in ti.tags, MANUAL_TEST_TAGS)
end
