# Multisegment well plotting functionality
# Requires GLMakie for interactive plotting

# NOTE: This module was developed with the assistance of AI agent tools.

"""
    compute_meters_drilled(well_model)

Compute meters drilled for each cell in a multisegment well.

# Arguments
- `well_model`: The multisegment well model (from a MultiModel)

# Returns
- Vector of meters drilled from the well start to each cell

Uses the neighborship matrix and cell centroids to compute cumulative distances
along the well path. Wells are treated as graphs where sections can connect
to any other section.
"""
function compute_meters_drilled(well_model)
    # Extract well data
    well_representation = well_model.domain.representation
    centroids = well_model.data_domain[:cell_centroids] # 3×n matrix
    neighborship = well_representation.neighborship # 2×m matrix [from_nodes; to_nodes]
    
    n_cells = size(centroids, 2)
    
    # Build graph adjacency list with distances
    graph = build_well_graph(centroids, neighborship)
    
    # Root/start node of the well (assumed to be the node with no incoming edges)
    start_node = 1
    
    # Compute shortest path distances from start to all nodes
    distances = compute_shortest_distances(graph, start_node, n_cells)
    
    return distances
end

"""
    build_well_graph(centroids, neighborship)

Build adjacency list representation of the well with edge weights as distances.
"""
function build_well_graph(centroids, neighborship)
    n_nodes = size(centroids, 2)
    n_segments = size(neighborship, 2)
    
    # Initialize adjacency list: node_id => [(neighbor_id, distance), ...]
    adj_list = [Vector{Tuple{Int, Float64}}() for _ in 1:n_nodes]
    
    # Add edges from neighborship matrix
    for i in 1:n_segments
        from_node = neighborship[1, i]
        to_node = neighborship[2, i]
        
        # Calculate Euclidean distance between centroids
        dist = euclidean_distance(centroids[:, from_node], centroids[:, to_node])
        
        # Add edge in both directions (undirected graph)
        push!(adj_list[from_node], (to_node, dist))
        # push!(adj_list[to_node], (from_node, dist))
    end
    
    return adj_list
end

"""
    euclidean_distance(p1, p2)

Compute Euclidean distance between two 3D points.
"""
function euclidean_distance(p1, p2)
    return norm(p1 - p2, 2)
end

"""
    find_well_start_node(neighborship, n_nodes)

Find the start/root node of the well.
"""
function find_well_start_node(neighborship, n_nodes)
    from_nodes = unique(neighborship[1, :])
    to_nodes = unique(neighborship[2, :])
    
    # Try to find nodes that only appear as 'from' (potential root nodes)
    potential_roots = setdiff(from_nodes, to_nodes)
    
    if !isempty(potential_roots)
        # Return the node with smallest index (conventional choice)
        return minimum(potential_roots)
    else
        # No clear root found (could be a cycle), default to node 1
        return 1
    end
end

"""
    compute_shortest_distances(graph, start_node, n_nodes)

Compute shortest path distances from start_node to all other nodes using Dijkstra's algorithm.
"""
function compute_shortest_distances(graph, start_node, n_nodes)
    # Initialize distances and visited status
    distances = fill(Inf, n_nodes)
    distances[start_node] = 0.0
    visited = fill(false, n_nodes)
    
    # Priority queue: vector of (distance, node_id) tuples
    pq = [(0.0, start_node)]
    
    while !isempty(pq)
        # Find and remove node with minimum distance
        min_idx = argmin([item[1] for item in pq])
        current_dist, current_node = splice!(pq, min_idx)
        
        if visited[current_node]
            continue
        end
        
        visited[current_node] = true
        
        # Update distances to all neighbors
        for (neighbor, edge_weight) in graph[current_node]
            if !visited[neighbor]
                new_dist = current_dist + edge_weight
                
                if new_dist < distances[neighbor]
                    distances[neighbor] = new_dist
                    # Add to priority queue
                    push!(pq, (new_dist, neighbor))
                end
            end
        end
    end
    
    # Replace any remaining Inf values with 0.0 (disconnected nodes)
    replace!(distances, Inf => 0.0)
    
    return distances
end

"""
    state_variable_unit_information(field::Symbol)

Return unit information for a state variable, similar to `well_target_information`.
Returns `missing` if the variable is not recognized.
"""
function state_variable_unit_information(field::Symbol)
    # Try well_target_information directly
    info = JutulDarcy.well_target_information(field)
    if !ismissing(info)
        return info
    end
    # Map common state variable names to unit information
    if field === :Pressure
        return JutulDarcy.well_target_information(Val(:bhp))
    elseif field === :Temperature
        return JutulDarcy.well_target_information(Val(:temperature))
    elseif field === :TotalMassFlux
        return JutulDarcy.well_target_information(
            symbol = :total_mass_flux,
            description = "Total mass flux",
            unit_type = :mass,
            unit_label = "kg/s",
            is_rate = true
        )
    else
        return missing
    end
end

"""
    convert_state_field_values(values, field::Symbol, unit_sys::String)

Convert state field values from SI to the specified unit system.
Returns `(converted_values, unit_label)`.
"""
function convert_state_field_values(values, field::Symbol, unit_sys::String)
    info = state_variable_unit_information(field)
    if ismissing(info)
        return values, string(field)
    end
    lbl = info.unit_label
    is_rate = info.is_rate
    factor = 1.0
    f_u, lbl = well_unit_conversion(unit_sys, lbl, info)
    if f_u isa Symbol
        values = convert_from_si.(values, f_u)
    else
        factor *= f_u
    end
    if is_rate
        lbl *= "/day"
        factor *= si_unit(:day)
    end
    values = factor .* values
    return values, lbl
end

"""
    plot_well_vs_meters_drilled!(ax, well_model, values; label=nothing, kwargs...)

Plot well data vs meters drilled on the given axis.

# Arguments
- `ax`: Makie axis to plot on
- `well_model`: The multisegment well model
- `values`: Vector of values (one per well cell) to plot

# Keyword arguments
- `label`: Base label for the plot
- `kwargs...`: Additional plotting arguments passed to lines!()

If the well has sections (`:section` key), each section is plotted with a distinct color.
Otherwise, plots all values as a single line.
"""
function JutulDarcy.plot_well_vs_meters_drilled!(ax, well_model, values;
    label=nothing, colors = missing, meters_drilled=missing,
    unit_sys="SI", field_name=missing, section_labels=missing, kwargs...)
    # Compute meters drilled for each cell if not provided
    if ismissing(meters_drilled)
        meters_drilled = compute_meters_drilled(well_model)
    end
    
    # Check if well has section data
    has_sections = haskey(well_model.data_domain, :section)

    n_cells = number_of_cells(well_model.data_domain)
    n_segments = number_of_faces(well_model.data_domain)
    if length(values) == n_cells
        # OK
    elseif length(values) == n_segments
        N = get_neighborship(well_model.data_domain.representation)
        meters_drilled = vec(sum(meters_drilled[N], dims=1)./2) # Average meters drilled for each segment
    else
        error("Length of values must match number of cells ($n_cells) or segments ($n_segments)")
    end

    # Apply unit conversion if field_name is provided
    if !ismissing(field_name)
        values, unit_label = convert_state_field_values(values, Symbol(field_name), unit_sys)
        ax.ylabel = unit_label
    end

    if has_sections
        sections = well_model.data_domain[:section]
        if length(values) == n_segments
            sections = max.(sections[N[1,:]], sections[N[2,:]])
        end
        num_sections = maximum(sections)
        if ismissing(colors)
            if num_sections <= 7
                colors = Makie.wong_colors(num_sections)  # Use Wong's color palette for better distinction
            elseif num_sections <= 10
                colors = cgrad(:tab10, num_sections, categorical = true) # Use Tab10 for up to 10 sections
            elseif num_sections <= 20
                colors = cgrad(:tab20, num_sections, categorical = true) # Use Tab20 for up to 20 sections
            else
                colors = cgrad(:BrBg, num_sections, categorical = true) # Use a diverging colormap for more than 20 sections (not ideal but better than repeating colors)
            end
        end
        
        # Plot each section with different colors
        for s in 1:num_sections
            cell_mask = sections .== s
            if any(cell_mask)  # Only plot if section has cells
                base_label = if !ismissing(section_labels)
                    section_labels[s]
                else
                    "Section $s"
                end
                section_label = if isnothing(label)
                    base_label
                else
                    "$label - $base_label"
                end
                
                lines!(ax, meters_drilled[cell_mask], values[cell_mask]; 
                      color=colors[s], label=section_label, linewidth=3, kwargs...)
            end
        end
    else
        # Plot single line without section coloring
        lines!(ax, meters_drilled, values; label=label, linewidth=3, kwargs...)
    end
    
    # Set axis labels
    ax.xlabel = "Meters Drilled"
    return ax
end

"""
    plot_well_vs_meters_drilled(well_model, values; figure_kwargs=NamedTuple(), axis_kwargs=NamedTuple(), kwargs...)

Create a new figure and plot well data vs meters drilled.

# Arguments
- `well_model`: The multisegment well model
- `values`: Vector of values (one per well cell) to plot

# Returns
- Figure handle

# Example
```julia
# Plot temperature vs meters drilled
fig = plot_well_vs_meters_drilled(well_model, temperature_values, 
                                  label="Temperature", 
                                  axis_kwargs=(ylabel="Temperature (°C)",))
```
"""
function JutulDarcy.plot_well_vs_meters_drilled(well_model, values; 
                                    figure_kwargs=NamedTuple(), 
                                    axis_kwargs=NamedTuple(), 
                                    unit_sys="SI",
                                    field_name=missing,
                                    section_labels=missing,
                                    kwargs...)
    fig = Figure(; figure_kwargs...)
    ax = Axis(fig[1, 1]; axis_kwargs...)
    JutulDarcy.plot_well_vs_meters_drilled!(ax, well_model, values; unit_sys=unit_sys, field_name=field_name, section_labels=section_labels, kwargs...)
    if haskey(well_model.data_domain, :section)
        Legend(fig[1, 2], ax, "Sections")
    end
    return fig
end

"""
    plot_well_states_interactive(well_model, states; 
                                  time=missing, names=missing, 
                                  resolution=(1200, 800), kwargs...)

Interactive GUI for plotting multisegment well data vs meters drilled.

# Arguments
- `well_model`: The multisegment well model
- `states`: Vector of states (one per timestep) or single state

# Keyword arguments
- `time`: Time vector corresponding to states
- `names`: Names for the different field variables
- `resolution`: Figure resolution
- `kwargs...`: Additional plotting arguments

# Features
- Field selection menu to choose which variable to plot
- Time stepping with slider for navigating through states
- Play button for movie-like animation
- Section-based coloring if well has sections
- Interactive controls for customizing the plot
"""
function JutulDarcy.plot_well_states_interactive(well_model, states; 
                                     time=missing, 
                                     names=missing, 
                                     unit_sys="Metric",
                                     section_labels=missing,
                                     resolution=(1200, 800), 
                                     kwargs...)
    
    # Handle single state case
    if states isa AbstractDict
        states = [states]
    end
    
    # Extract available fields from first state
    first_state = states[1]
    field_names = collect(keys(first_state))
    
    # Filter to numeric fields with correct dimensions
    meters_drilled = compute_meters_drilled(well_model)
    n_cells = length(meters_drilled)
    n_segments = number_of_faces(well_model.data_domain)
    
    valid_fields = String[]
    for field in field_names
        val = first_state[field]
        if val isa AbstractVector && eltype(val) <: Number
            if length(val) ∈ [n_cells, n_segments]
                push!(valid_fields, string(field))
            end
        end
    end
    
    if isempty(valid_fields)
        error("No valid fields found for plotting")
    end
    
    # Compute global y-limits for consistent scaling across timesteps
    field_y_limits = Dict{String, Tuple{Float64, Float64}}()
    for field in valid_fields
        field_symbol = Symbol(field)
        all_values = Float64[]
        for state in states
            append!(all_values, state[field_symbol])
        end
        y_min, y_max = extrema(all_values)
        # Add 5% padding
        y_range = y_max - y_min
        y_padding = max(y_range * 0.05, eps(Float64))
        field_y_limits[field] = (y_min - y_padding, y_max + y_padding)
    end
    
    # Setup observables for reactive GUI
    n_states = length(states)
    state_index = Observable{Int}(1)
    field_selection = Observable{String}(valid_fields[1])
    unit_sys_obs = Observable{String}(unit_sys)
    is_playing = Observable{Bool}(false)
    play_speed = Observable{Float64}(1.0)  # States per second
    
    # Create figure and layout
    fig = Figure(size=resolution)
    
    # Main plot axis
    ax = Axis(fig[1:4, 1:3], 
              xlabel="Meters Drilled", 
              ylabel="Field Value",
              aspect=AxisAspect(2))
    
    # Control panel at bottom
    control_layout = fig[5, :] = GridLayout()
    
    # Field selection menu
    field_menu = Menu(fig, options=valid_fields, default=valid_fields[1], width=150)
    on(field_menu.selection) do s
        field_selection[] = s
    end

    # Unit system menu
    unit_menu = Menu(fig, options=["Metric", "SI", "Field"], default=unit_sys_obs[], width=100)
    on(unit_menu.selection) do s
        unit_sys_obs[] = s
    end
    
    # Time stepping controls
    time_slider = Slider(fig, range=1:n_states, value=1, snap=true)
    
    # Manual bidirectional connection to avoid circular dependency
    # Flag to prevent infinite loops during programmatic updates
    updating_slider = Ref(false)
    
    # Slider → state_index
    on(time_slider.value) do val
        if !updating_slider[]
            state_index[] = val
        end
    end
    
    # state_index → Slider (for programmatic updates like play button)
    on(state_index) do val
        updating_slider[] = true
        try
            set_close_to!(time_slider, val)
        finally
            updating_slider[] = false
        end
    end
    
    # Play/pause button
    play_button = Button(fig, label="Play", width=80)
    
    # Speed control
    speed_input = Textbox(fig, placeholder="1.0", width=60, validator=s -> tryparse(Float64, s) !== nothing)
    speed_input.stored_string = "1.0"
    
    on(speed_input.stored_string) do s
        val = tryparse(Float64, s)
        if val !== nothing && val > 0
            play_speed[] = val
        end
    end
    
    # Reset button
    reset_button = Button(fig, label="Reset", width=60)
    on(reset_button.clicks) do _
        state_index[] = 1
        is_playing[] = false
        play_button.label = "Play"
    end
    
    # Axis reset buttons
    reset_x_button = Button(fig, label="Reset X", width=70)
    reset_y_button = Button(fig, label="Reset Y", width=70)
    
    on(reset_x_button.clicks) do _
        # Reset x-axis to full range of meters drilled
        x_min, x_max = extrema(meters_drilled)
        x_range = x_max - x_min
        x_padding = max(x_range * 0.02, eps(Float64))
        xlims!(ax, x_min - x_padding, x_max + x_padding)
    end
    
    on(reset_y_button.clicks) do _
        # Reset y-axis to full range of current field (converted to current unit system)
        field = field_selection[]
        y_min, y_max = field_y_limits[field]
        lims, _ = convert_state_field_values([y_min, y_max], Symbol(field), unit_sys_obs[])
        ylims!(ax, lims[1], lims[2])
    end
    
    # Layout controls
    control_layout[1, 1] = Label(fig, "Field:", halign=:right)
    control_layout[1, 2] = field_menu
    control_layout[1, 3] = Label(fig, "Units:")
    control_layout[1, 4] = unit_menu
    control_layout[1, 5] = Label(fig, "Time step:")
    control_layout[1, 6] = time_slider
    control_layout[1, 7] = play_button
    control_layout[1, 8] = Label(fig, "Speed:")
    control_layout[1, 9] = speed_input
    control_layout[1, 10] = reset_button
    control_layout[1, 11] = reset_x_button
    control_layout[1, 12] = reset_y_button
    
    # Animation timer
    animation_timer = Timer(0.1)  # Will be started/stopped as needed
    
    # Play button functionality
    on(play_button.clicks) do _
        if is_playing[]
            # Stop playing
            is_playing[] = false
            play_button.label = "Play"
            close(animation_timer)
        else
            # Start playing
            is_playing[] = true 
            play_button.label = "Pause"
            
            # Create new timer for animation
            dt = 1.0 / play_speed[]
            animation_timer = Timer(dt, interval=dt)
            
            @async while is_playing[] && isopen(animation_timer)
                wait(animation_timer)
                if state_index[] < n_states
                    state_index[] = state_index[] + 1
                else
                    # End of animation
                    is_playing[] = false
                    play_button.label = "Play"
                    close(animation_timer)
                    break
                end
            end
        end
    end
    
    # Update speed during playback
    on(play_speed) do speed
        if is_playing[]
            close(animation_timer)
            dt = 1.0 / speed
            animation_timer = Timer(dt, interval=dt)
        end
    end
    
    # Main plotting logic - reactive to state and field changes
    function update_plot()
        idx = state_index[]
        field = field_selection[]
        usys = unit_sys_obs[]
        
        # Clear existing plots
        empty!(ax)
        
        # Get current state and field data  
        current_state = states[idx]
        field_data = Float64.(current_state[Symbol(field)])
        
        # Convert units
        field_data, unit_label = convert_state_field_values(field_data, Symbol(field), usys)
        
        # Plot data vs meters drilled (with section support)
        JutulDarcy.plot_well_vs_meters_drilled!(ax, well_model, field_data; meters_drilled=meters_drilled, section_labels=section_labels, kwargs...)
        
        # Set consistent y-axis limits (convert from SI)
        y_min, y_max = field_y_limits[field]
        lims, _ = convert_state_field_values([y_min, y_max], Symbol(field), usys)
        ylims!(ax, lims[1], lims[2])
        
        # Update title
        title_str = if n_states > 1
            if ismissing(time)
                "$field - Step $idx/$n_states"
            else
                if time[idx] isa DateTime
                    "$field - $(time[idx])"
                else
                    "$field - Time: $(round(time[idx], digits=2))"
                end
            end
        else
            field
        end
        ax.title = title_str
        
        # Update axis label with unit
        ax.ylabel = unit_label
    end
    
    # Listen to state, field, and unit system observables
    on(state_index) do _
        update_plot()
    end
    
    on(field_selection) do _
        update_plot()
    end
    
    on(unit_sys_obs) do _
        update_plot()
    end
    
    # Information panel (top right)
    info_layout = fig[1, 4] = GridLayout()
    
    # State info
    state_label = @lift("State: $($state_index)/$n_states")
    info_layout[1, 1] = Label(fig, state_label)
    
    if !ismissing(time) && length(time) == n_states
        time_label = if time[1] isa DateTime
            @lift("Time: $(time[$state_index])")
        else
            @lift("Time: $(round(time[$state_index], digits=3))")
        end
        info_layout[2, 1] = Label(fig, time_label)
    end
    
    # Add section legend to the right of the plot
    has_sections = haskey(well_model.data_domain, :section)
    if has_sections
        n_sections = maximum(well_model.data_domain[:section])
        # Determine section colors (matching logic in plot_well_vs_meters_drilled!)
        kw = Dict{Symbol,Any}(kwargs)
        user_colors = get(kw, :colors, missing)
        if !ismissing(user_colors)
            section_colors = user_colors
        elseif n_sections <= 7
            section_colors = Makie.wong_colors(n_sections)
        elseif n_sections <= 10
            section_colors = cgrad(:tab10, n_sections, categorical = true)
        elseif n_sections <= 20
            section_colors = cgrad(:tab20, n_sections, categorical = true)
        else
            section_colors = cgrad(:BrBg, n_sections, categorical = true)
        end
        entries = [LineElement(color = section_colors[s], linewidth = 3) for s in 1:n_sections]
        labels = if !ismissing(section_labels)
            [section_labels[s] for s in 1:n_sections]
        else
            ["Section $s" for s in 1:n_sections]
        end
        Legend(fig[2:4, 4], entries, labels, "Sections")
    end
    
    # Trigger initial plot
    notify(state_index)
    
    return fig
end

"""
    plot_well_states_interactive(well_name, model, states; kwargs...)

Convenience function to plot well states by extracting well model from a MultiModel.

# Arguments
- `well_name`: Symbol name of the well (e.g., :Well1)
- `model`: MultiModel containing the well
- `states`: Vector of states or simulation results

# Example
```julia
# Plot all timesteps interactively
fig = plot_well_states_interactive(:AGS_supply, case.model, results.result.states)

# Plot specific states with time information
times = convert_from_si.(cumsum(case.dt), :day)
fig = plot_well_states_interactive(:AGS_supply, case.model, results.result.states[1:10],
                                  time=times[1:10])
```
"""
function JutulDarcy.plot_well_states_interactive(well_name::Symbol, model, states; kwargs...)
    well_model = model.models[well_name]
    states_for_well = [state[well_name] for state in states]
    return JutulDarcy.plot_well_states_interactive(well_model, states_for_well; kwargs...)
end