# coding=utf-8

import numpy as np


#############   plane equasion AND PROJECTION COORDINATES ################
def plane_eq(x, y, z):
    """This function will return the plane equasion parameters a,b,c and d
        for each three dots that are in the plane """

    # defining two vectors that start from the same dot, let's choose dot x
    p_xy = y - x
    p_xz = z - x
    # finding ortogonal vector to the ones we have just created
    norm_p = np.cross(p_xy, p_xz)
    # for any dot in the plane let's choose X, norm_p * X = d
    d = -1 * np.dot(norm_p, x)
    # direction vector presents a,b and c
    a = norm_p[0]
    b = norm_p[1]
    c = norm_p[2]

    return a, b, c, d


#############################################################################

########### PROJECTION IN PLANE AND INTERACTING RING ATOM COORDINATES #############
def projection_in_plane(A, B, C, xs, ys, zs):
    """ This function takes three dots A, B, C to form a plane and
    returns coordinates of projection of atom_O in the plane
    A,B,C = arrays
    atom_O = list with one array

    For further purposes this function returns coordinates of oxigen projection (coord_proj).
    """
    (a, b, c, d) = plane_eq(A, B, C)

    # parameter t for parametric equasion is
    part1_t = a * xs + b * ys + c * zs + d
    part2_t = a * a + b * b + c * c
    part3_t = part1_t / part2_t
    t = -1 * part3_t

    x_proj = a * t + xs
    y_proj = b * t + ys
    z_proj = c * t + zs

    coord_proj = []
    coord_proj.append(x_proj)
    coord_proj.append(y_proj)
    coord_proj.append(z_proj)

    return coord_proj


##########################################################################



if __name__ == '__main__':

    xd = """
ATOM  12407  CHD HEM     2      -5.164 -30.729  20.121  1.00 21.79      EHEM
    """

    xp = np.array([-1.505, -33.622,  21.512])
    yp = np.array([-3.669, -34.209,  25.866])
    zp = np.array([-7.350, -31.365,  24.419])

    xfe= -4.520
    yfe= -32.644
    zfe = 22.892

    fe = np.array([xfe, yfe, zfe])

    point = projection_in_plane(xp, yp, zp, xfe, yfe, zfe)

    print point

    vector = np.subtract(np.array(point),fe)

    vector_mag = np.linalg.norm(vector)
    unitvector = vector / vector_mag

    # C_atom = point + unitvector * 1.7 #position above the prjection for C
    # O_atom = point + unitvector * 2.83 #position above the prjection for O

    C_atom = fe + unitvector * 1.7 #position above the prjection for C
    O_atom = fe + unitvector * 2.83 #position above the prjection for O

    print C_atom
    print O_atom



    


















