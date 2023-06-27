# Generate test geometry but just its base (i.e. just the circle) specifying h# Generation of test geometry and producing its triangular mesh
# I also generated a triangular mesh in a square because I'm going to 
# first test my fmm method here (very naivly)
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import pi, sqrt, floor, ceil
import meshpy.triangle as triangle
from matplotlib.patches import Arc
import numpy.linalg as la
from matplotlib.colors import ListedColormap
import colorcet as cc



colormap2 = "cet_linear_worb_100_25_c53_r"
colormap2_r = "cet_linear_worb_100_25_c53"




plt.ion()

h =  0.01
h_string = str(h)
eta1 = 1.0
eta2 = 1.45
currentH = "H9"
path_figures = "/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/"
path_info = '/Users/marianamartinez/Documents/Curvy-JMM/JMM/'

############# Generate the two circles

my_dpi=96

def circ_base(t):
    '''
    Parametrization of the base cirlce (base of the snowman) for the test geometry.
    t can be a number or an array or numbers. If t is a number the output are the
    coordinates [x, y]. If t is an array the output is a n, 2 array, the first 
    column corresponds to the coordinates of x, the second to the coordinates of y.
    '''
    points = np.zeros((len(t), 2))
    tan_points = np.zeros((len(t), 2))
    points[:, 0] = 10*np.cos(t)
    points[:, 1] = 10*np.sin(t)
    tan_points[:, 0] = -np.sin(t)
    tan_points[:, 1] = np.cos(t)
    return points, tan_points



############# Useful functions

def round_trip_connect(start, end, boundary_tan = False, tans = None, hc = None):
    facets = [(i, i + 1) for i in range(start, end)] + [(end, start)]
    if( boundary_tan ):
        # This means that we also want to get the boundary tan struct
        B = []
        for i in range(start, end):
            tans_i = tans[i]/norm(tans[i])*hc[i]
            tans_i1 = tans[i+1]/norm(tans[i+1])*hc[i+1]
            B.append( (i, i+1, tans_i[0], tans_i[1], tans_i1[0], tans_i1[1] ) )
        tans_n = tans[end]/norm(tans[end])*hc[end]
        tans_1 = tans[start]/norm(tans[start])*hc[start]
        B.append( (start, end, tans_1[0], tans_1[1], tans_n[0], tans_n[1] ) )
        return facets, B
    else:
        return facets

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]

def needs_refinement(vertices, area):
    max_area = h
    return bool(area > max_area)
    
def nPoints(h):
    '''
    With this function we know how many points we need in the top of the snowman and in the bottom, depending on which h we want
    WE ASSUME THAT WE ARE DISTRIBUTING UNIFORMLY THESE POINTS
    '''
    nBottom = ceil(15*pi/h) 
    return nBottom

Nbase = nPoints(h)

t = np.linspace(3/4, 11/4, num = Nbase , endpoint = False)
t = np.multiply(t, pi)

points_base, tan_base = circ_base(t)
hc = [norm( points_base[0, :] - points_base[1, :] )]*len(points_base) # All are the same because the points are equispaced on the circle
    

# The points on the outline of the snowman and the tangent
points =  [(-15, -10)] # start with x0 so that in each mesh we know its index
hc = [0] + hc
points.extend([ (points_base[i, 0], points_base[i, 1]) for i in range( len(points_base) ) ] )
# Add the facets
tans = np.zeros((len(points), 2))
tans[1:, :] = tan_base
facets, mesh_square_boundaryCurve = round_trip_connect(1, len(points) -1, boundary_tan = True, tans = tans, hc = hc ) # we don't want to connnect x0
mesh_square_boundaryCurve = np.array(mesh_square_boundaryCurve)


##############################################################################
##############################################################################
##############################################################################
##############################################################################
###### SQUARE ON THE OUTSIDE

edges_square = [(-18, -18), (18, -18), (18, 24), (-18, 24)]

# Add them to the points 
points_square = points + edges_square
# Add the facets
facets_square = facets + round_trip_connect(len(points), len(points_square)-1) 

# Add them to the triangulation (this is a new triangulation)
info_square = triangle.MeshInfo()
info_square.set_points(points_square)
info_square.set_facets(facets_square)

mesh_square = triangle.build(info_square,volume_constraints= True, refinement_func=needs_refinement, attributes=True, generate_faces=True)
mesh_square_points = np.array(mesh_square.points) # These are the points we want to export
mesh_square_tris = np.array(mesh_square.elements) # These are thee faces we want to export
mesh_square_tris = np.sort(mesh_square_tris, axis = 1) # But we want them in order
mesh_square_neigTriangles = np.array(mesh_square.neighbors) # These are the neighbors we want to export
mesh_square_neigTriangles = np.sort(mesh_square_neigTriangles, axis = 1) # We also want these sorted
mesh_square_facets = np.array(mesh_square.facets) # These are the edges we want to exportmesh_square_boundaryCurve = np.array(mesh_square_boundaryCurve) # Because our C code adds the NULL automatically
mesh_square_facets = np.sort(mesh_square_facets, axis = 1) # We also want these sorted

print("Triangulation done, number of points: ", len(mesh_square_points))

# Now we want to create an np array where the rows are the #of the point, the columns each one of its neighbours
# This is the most naive way of creating such thing, might be useful to optimize it later (?)
# We look in the row of the mesh.elements file
N_points_square = len(mesh_square_points)
mesh_square_neigh = [[] for i in range(N_points_square)]
MaxN_square = 0


# Create the list of lists

for t in range(len(mesh_square_tris)):
    # We iterate through each triangle
    p0 = mesh_square_tris[t,0]
    p1 = mesh_square_tris[t,1]
    p2 = mesh_square_tris[t,2]
    if p0 not in mesh_square_neigh[p1]:
        mesh_square_neigh[p1].append(p0)
    if p0 not in mesh_square_neigh[p2]:
        mesh_square_neigh[p2].append(p0)
    if p1 not in mesh_square_neigh[p0]:
        mesh_square_neigh[p0].append(p1)
    if p1 not in mesh_square_neigh[p2]:
        mesh_square_neigh[p2].append(p1)
    if p2 not in mesh_square_neigh[p0]:
        mesh_square_neigh[p0].append(p2)
    if p2 not in mesh_square_neigh[p1]:
        mesh_square_neigh[p1].append(p2)

print("List of neighbors done")


# Create the list for the incident faces and edges

mesh_IncidentFaces_sq = [[] for i in range(N_points_square)]
mesh_square_facets = [list(f) for f in mesh_square_facets]

for t in range(len(mesh_square_tris)):
    # Iterate through each triangle
    p0 = mesh_square_tris[t, 0]
    p1 = mesh_square_tris[t, 1]
    p2 = mesh_square_tris[t, 2]
    # For the incident faces
    mesh_IncidentFaces_sq[p0].append(t)
    mesh_IncidentFaces_sq[p1].append(t)
    mesh_IncidentFaces_sq[p2].append(t)
    # For the edges, all of of them, then we will filter them
    mesh_square_facets.append([p0, p1])
    mesh_square_facets.append([p0, p2])
    mesh_square_facets.append([p1, p2])

# Format the way we need it
mesh_square_facets = np.array(mesh_square_facets)
mesh_square_facets = np.sort(mesh_square_facets, axis = 1)
mesh_square_facets = np.unique(mesh_square_facets, axis = 0)



print("List of incident faces done")
print("List of edges done")



# Now we create the np array for the edges in face
edgesInFaces = np.ones((len(mesh_square_tris), 4), dtype = np.int64 )

for i in range(len(mesh_square_facets)):
    # Indices
    indFrom = np.where(mesh_square_tris == mesh_square_facets[i, 0])[0]
    indTo = np.where(mesh_square_tris == mesh_square_facets[i, 1])[0]
    indFaces = list(set(indFrom) & set(indTo) ) # This should be a list of length 2 or 1
    assert( len(indFaces) <= 2)
    for k in range(len(indFaces)):
        j = edgesInFaces[indFaces[k], 0]
        edgesInFaces[indFaces[k], j] = i
        edgesInFaces[indFaces[k], 0] += 1

edgesInFaces = edgesInFaces[:, 1:]

print("List edges in face done")

# We need to put mesh_square_boundaryCurve in order
for i in range(len(mesh_square_boundaryCurve)):
        # Indices
        indFrom = np.where(mesh_square_facets == mesh_square_boundaryCurve[i, 0])[0]
        indTo = np.where(mesh_square_facets == mesh_square_boundaryCurve[i, 1])[0]
        indEdge = list(set(indFrom) & set(indTo) ) # This should be a list of length 1
        assert( len(indEdge) == 1)
        mesh_square_boundaryCurve[i, 1] = indEdge[0]


mesh_square_boundaryCurve = mesh_square_boundaryCurve[:, 1:]
print(" Mesh square boundary organized")

    
# Create the labels per face


faces_label = []
colors = []
r = 5*sqrt(2)
for fi in range(len(mesh_square_tris)):
    # we need to get the coordinates of the 3 points that define that upper triangle
    p1 =  mesh_square_points[mesh_square_tris[fi, 0] , :]
    p2 =  mesh_square_points[mesh_square_tris[fi, 1] , :]
    p3 =  mesh_square_points[mesh_square_tris[fi, 2] , :]
    if ( sqrt(p1[0]**2 + p1[1]**2)<= 10 ) and ( sqrt(p2[0]**2 + p2[1]**2)<= 10 ) and (sqrt(p3[0]**2 + p3[1]**2)<= 10):
        faces_label += [eta2]
        colors += [0.3]
    elif( (sqrt(p1[0]**2 + p1[1]**2) + sqrt(p2[0]**2 + p2[1]**2) + sqrt(p3[0]**2 + p3[1]**2))<= 30  ):
        faces_label += [eta2]
        colors += [0.3]
    else:
        faces_label += [eta1]
        colors += [0.1]

print("Lists of faces labels done")




# #Plot
plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.title('Delaunay triangulation of test geometry, H = '+h_string)
plt.savefig(path_figures + currentH + '/'+ currentH + '_TriangulationWhite.png', dpi=my_dpi * 10)
plt.xlim([-18, 18])
plt.ylim([-18, 24])
plt.show(block=False)


viridis = plt.get_cmap(colormap2, 256)
new_colors = viridis(np.linspace(0, 1, 256))
blue = new_colors[0]
purple = new_colors[200]
mid = floor(len(new_colors)/2)
new_colors = [blue]*mid
new_colors += [purple]*mid
newcmp = ListedColormap(new_colors)


plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
plt.gca().set_aspect('equal')
plt.tripcolor(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, faces_label, cmap = colormap2)
plt.triplot(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, '-.', lw=0.5, c='#00fffb')
plt.title('Delaunay triangulation of test geometry, H = '+h_string)
plt.xlim([-18, 18])
plt.ylim([-18, 24])
plt.show(block = False)
plt.savefig(path_figures + currentH + '/'+ currentH + '_Triangulation.png', dpi=my_dpi * 10)


# With the tangents
plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
plt.gca().set_aspect('equal')
plt.tripcolor(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, faces_label, cmap = colormap2)
plt.triplot(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, '-.', lw=0.5, c='#00fffb')
plt.quiver(points_base[:,0], points_base[:,1], mesh_square_boundaryCurve[:, 1],
           mesh_square_boundaryCurve[:, 2], width = 0.005, scale = 1.0, scale_units = 'xy')
plt.title('Delaunay triangulation of test geometry, H = '+h_string)
plt.xlim([-18, 18])
plt.ylim([-18, 24])
plt.show(block = False)
plt.savefig(path_figures + currentH + '/'+ currentH + '_TriangulationTangents.png', dpi=my_dpi * 10)


## Save them as well

np.savetxt(path_info + currentH + '/'+ currentH +'_Edges.txt', mesh_square_facets.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_MeshPoints.txt', mesh_square_points, delimiter =', ', fmt = '%.8f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_Faces.txt', mesh_square_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_NeighTriangles.txt', mesh_square_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')


# Sort mesh_square_boundaryCurve
IndOrg = np.argsort( mesh_square_boundaryCurve[:,0] ) # Order according to zeroth column (i.e. edge number)
mesh_square_boundaryCurve = mesh_square_boundaryCurve[IndOrg]
np.savetxt(path_info + currentH + '/'+ currentH +'_BoundaryCurve.txt', mesh_square_boundaryCurve, delimiter =', ', fmt = '%.8f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_Indices.txt', np.array(faces_label), delimiter =', ', fmt = '%.8f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_EdgesInFace.txt', edgesInFaces.astype(int), delimiter =', ', fmt ='%.0f' )

separator = "," 

with open(path_info + currentH + '/'+ currentH + "_Neigh.txt", "w") as out_file:
    for l in mesh_square_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open(path_info + currentH + '/'+ currentH + "_IncidentFaces.txt", "w") as out_file:
    for l in mesh_IncidentFaces_sq:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
plt.show()
