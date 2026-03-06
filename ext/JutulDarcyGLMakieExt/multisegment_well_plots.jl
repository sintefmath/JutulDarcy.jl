# Multisegment well plotting functionality
# Requires GLMakie for interactive plotting

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
    
    # Find the root/start node of the well
    start_node = find_well_start_node(neighborship, n_cells)
    
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
        push!(adj_list[to_node], (from_node, dist))
    end
    
    return adj_list
end

"""
    euclidean_distance(p1, p2)

Compute Euclidean distance between two 3D points.
"""
function euclidean_distance(p1, p2)
    return norm(p1 - p2)
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
function JutulDarcy.plot_well_vs_meters_drilled!(ax, well_model, values; label=nothing, kwargs...)
    # Compute meters drilled for each cell
    meters_drilled = compute_meters_drilled(well_model)
    
    # Check if well has section data
    has_sections = haskey(well_model.data_domain, :section)
    
    if has_sections
        sections = well_model.data_domain[:section]
        num_sections = maximum(sections)
        # colors = cgrad(:BrBg, num_sections, categorical = true)
        colors = Makie.wong_colors(num_sections)  # Use Wong's color palette for better distinction
        
        # Plot each section with different colors
        for s in 1:num_sections
            cell_mask = sections .== s
            if any(cell_mask)  # Only plot if section has cells
                section_label = if isnothing(label)
                    "Section $s"
                else
                    "$label - Section $s"
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
                                    kwargs...)
    fig = Figure(; figure_kwargs...)
    ax = Axis(fig[1, 1]; axis_kwargs...)
    JutulDarcy.plot_well_vs_meters_drilled!(ax, well_model, values; kwargs...)
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
    
    valid_fields = String[]
    for field in field_names
        val = first_state[field]
        if val isa AbstractVector && length(val) == n_cells && eltype(val) <: Number
            push!(valid_fields, string(field))
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
    is_playing = Observable{Bool}(false)
    play_speed = Observable{Float64}(1.0)  # States per second
    
    # Create figure and layout
    fig = Figure(size=resolution)
    
    # Main plot axis
    ax = Axis(fig[1, 1:3], 
              xlabel="Meters Drilled", 
              ylabel="Field Value")
    
    # Control panel at bottom
    control_layout = fig[2, :] = GridLayout()
    
    # Field selection menu
    field_menu = Menu(fig, options=valid_fields, default=valid_fields[1], width=150)
    on(field_menu.selection) do s
        field_selection[] = s
        # Update axis label
        ax.ylabel = s
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
        # Reset y-axis to full range of current field
        field = field_selection[]
        y_min, y_max = field_y_limits[field]
        ylims!(ax, y_min, y_max)
    end
    
    # Layout controls
    control_layout[1, 1] = Label(fig, "Field:", halign=:right)
    control_layout[1, 2] = field_menu
    control_layout[1, 3] = Label(fig, "Time step:")
    control_layout[1, 4] = time_slider
    control_layout[1, 5] = play_button
    control_layout[1, 6] = Label(fig, "Speed:")
    control_layout[1, 7] = speed_input
    control_layout[1, 8] = reset_button
    control_layout[1, 9] = reset_x_button
    control_layout[1, 10] = reset_y_button
    
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
        
        # Clear existing plots
        empty!(ax)
        
        # Get current state and field data  
        current_state = states[idx]
        field_data = current_state[Symbol(field)]
        
        # Plot data vs meters drilled (with section support)
        JutulDarcy.plot_well_vs_meters_drilled!(ax, well_model, field_data)
        
        # Set consistent y-axis limits
        y_min, y_max = field_y_limits[field]
        ylims!(ax, y_min, y_max)
        
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
        
        # Update axis label
        ax.ylabel = field
    end
    
    # Listen to both observables
    on(state_index) do _
        update_plot()
    end
    
    on(field_selection) do _
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
    
    # Check for sections
    has_sections = haskey(well_model.data_domain, :section)
    if has_sections
        n_sections = maximum(well_model.data_domain[:section])
        info_layout[3, 1] = Label(fig, "Sections: $n_sections")
        info_layout[4, 1] = Label(fig, "Colors by section")
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