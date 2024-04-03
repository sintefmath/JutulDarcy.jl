cart_dims = (100, 1, 50)
physical_dims = (1000.0, 1.0, 100.0)
mesh = UnstructuredMesh(CartesianMesh(cart_dims, physical_dims))

points = mesh.node_points
for (i, pt) in enumerate(points)
    x, y, z = pt
    x_u = 2*π*x/1000.0
    dz = 0.05*x + 30*sin(2.0*x_u) + 20*sin(5.0*x_u) + 10*sin(10.0*x_u) + 5*sin(25.0*x_u)
    points[i] = pt + [0, 0, dz]
end

using GLMakie
plot_mesh(mesh)