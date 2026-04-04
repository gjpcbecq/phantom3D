"""
Tools code for the phantom3D 

Author: GJC Becq
Date: 2026-04-03
Copyright: CeCILL (compatible GNU see LICENSE.txt)

"""

import ext
    
def test_flips_orientations(flips1, flips2, permutations1, permutations2): 
    """
    s = test_flips_orientations(flips1, flips2, permutations1, permutations2) 

    compare flips1 and flips2, permutations1 and permutations2
    return a string corresponding to failed or passed tests between flips and orientations
    
    
    flips1: tuple, (flip1, flip2, flip3) with flip1 is in {-1, 1}
    permutations1: tuple, indeed a permutation example (0, 1, 2)

    but elements can be any tuples of 3 elements.
    
    input

    flips1: tuple, (flip1, flip2, flip3) with flip1 is in {-1, 1}
    flips2: tuple, same as flips1
    permutations1: tuple, indeed a permutation example (0, 1, 2)
    permutations2: tuple, same as permutations1

    output

    s: a string, containing the differences between flips1 and flips2, permutations1 and permutations2. A "F" is added to the string if differences are observed btween elemens of the tuples, a "." otherwise.  
    
    """
    
    s = ""
    for i in range(3): 
        if flips1[i] == flips2[i]: 
            s += "." 
        else: 
            s += "F"
    for i in range(3): 
        if permutations1[i] == permutations2[i]: 
            s += "." 
        else: 
            s += "F"
    return s
        

def compute_A_inside_B_criterion_with_translation(A, B, i, j, k, I1, J1, K1):
    """
    c = compute_A_inside_B_criterion_with_translation(A, B, i, j, k, I1, J1, K1)

    compute the criterion (cosine similarity score) of A inside B.
    A is flatten(), and extract of B, which is the size of A is also flatten() before computing the cosine similarity score. 
    The extract of B corresponds to the occurence of A at shift (i, j, k). 
    
    input
    
    A: ndarray with shape (I1, J1, K1)
    B: ndarray with shape (I2, J2, K2)
    I1, J1, K1: int, size of the volume A inside B so that I2 >= I1, J2 >= J1, K2 >= K1
    i, j, k: int, the indices of the shift (0 no shift)
    0 <= i <= I2 - I1
    0 <= j <= J2 - J1
    0 <= k <= K2 - K1

    output

    c: the similarity score of A inside B for the shift (i,j,k)
    
    """
    s1 = slice(i, i + I1)
    s2 = slice(j, j + J1)
    s3 = slice(k, k + K1)
    Bextract = B[s1, s2, s3]
    c = compute_criterion(Bextract.flatten(), A.flatten())
    return c

def compute_A_inside_B_criterion(A, B): 
    """
    C = compute_A_inside_B_criterion(A, B)

    A should be smaller than B. (I1 <= I2, J1 <= J2, K1 <= K2) 
    All shifts are evaluated to compute the score. 
    
    Input

    A: ndarray (I1, J1, K1), template
    B: ndarray (I2, J2, K2), analyzed volume

    Output

    C: ndarray (I, J, K), contains the cosine similarity score evaluate at each translations (i,j,k) $(i,j,k) \\in (I,J,K)$   
    
    """
    
    (I1, J1, K1) = A.shape
    (I2, J2, K2) = B.shape
    I = I2 - I1 + 1
    J = J2 - J1 + 1
    K = K2 - K1 + 1
    C = ext.np.zeros((I, J, K))
    for i in range(I): 
        for j in range(J): 
            for k in range(K): 
                C[i, j, k] = compute_A_inside_B_criterion_with_translation(A, B, i, j, k, I1, J1, K1)
    return C

def compute_A_inside_B_best_translation(A, B): 
    """
    (ix, jx, kx, cx) = compute_A_inside_B_best_translation(A, B)

    input

    A: ndarray (I1, J1, K1), template
    B: ndarray (I2, J2, K2), analyzed volume

    output

    (ix, jx, kx, cx) : tuple, ix, jx, kx are int and corresponds to the indices of the best shift. cx is a float $\\in [-1, 1]$ and corresponds to the best cosine similarity score obtained at ix, jx, kx.   

    see also
    C = compute_A_inside_B_criterion(A, B)
    
    """
    C = compute_A_inside_B_criterion(A, B)
    imax = ext.np.argmax(ext.np.abs(C))
    (ix, jx, kx) = ext.np.unravel_index(imax, C.shape)
    cx = C[ix, jx, kx]
    return (ix, jx, kx, cx)    

def compute_criterion(x, xref, eps=1e-5): 
    """
    c = compute_criterion(x, xref, eps=1e-5): 
    compute the criterion here the cosine similarity between vector x and xref. 

    $$ C = (x.T\\, xref) / (||x||\\, ||xref||)^{0.5}$$ 
    
    input 
    
    x: ndarray of shape (n,) 
    xref: ndarray of shape (n,)

    output
    
    c: float, criterion score, cosine similarity score (-1, 1)
    
    """
    nref = xref.dot(xref)
    n = x.dot(x)
    if n < eps: 
        return 0
    den = (nref * n) ** 0.5
    c = xref.dot(x) / den
    return c

def compute_criterion_A_Bx(A, Bx): 
    """
    c = compute_criterion_A_Bx(A, Bx)

    compute the criterion (cosine similarity score between volume A and Bx)
    Bx is an array extracted from B such that Bx is the same shape as A

    A and Bx are flatten 
    
    input

    A: ndarray of shape (I,J,K)
    Bx: ndarray of shape (I,J,K)
    
    output

    c: float (-1, 1) similarity score 
    
    """
    c = compute_criterion(Bx.flatten(), A.flatten())
    return c

def compute_criterion_corrcoef(x, xref, eps=1e-5): 
    """
    c = compute_criterion_corrcoef(x, xref, eps=1e-5): 

    compute the correlation coefficient (See Pearson correlation)

    
    input 
    
    x: ndarray of shape (n,) 
    xref: ndarray of shape (n,)

    output
    
    c: float, criterion score, corrcoef score (-1, 1)
    
    
    """
    nref = xref.dot(xref.T)
    n = x.dot(x.T)
    if n < eps: 
        return 0
    c = ext.np.corrcoef(xref, x)[0, 1]
    return c

def compute_flips_permutations(aff1):
    """
    (flips, permutations) = compute_flips_permutations(aff1)

    given aff1, the maximal value is detectd on each axis, the sign of this maximal value is also detected, then the correspondance with flips and permutations is determined. 

    input

    aff1: ndarray, an affine transform of shape (4, 4)
    
    output

    flips: tuple of boolean (flips for axes i, j and k) example (True, True, False)
    permutations: tuple of int in [0, 1, 2] corresponding to permutations of axes. 
        Example (2, 1, 0) is the permutation of axes 0 and 2, with axis 1 not permuted. 

    This is a fast check for axes aligned with few rotations, I don't know what it does when rotations are present in the matrix. 
    
    """
    flips = [False, False, False]
    permutations = [0, 1, 2]
    for i in range(3): 
        j = ext.np.argmax(ext.np.abs(aff1[i, :3]))
        permutations[i] = j
        if aff1[i, j] < 0: 
            flips[i] = True
        else: 
            flips[i] = False
    return (flips, permutations)

def compute_new_pattern(Ap0, flip, permutation): 
    """
    Ap1 = compute_new_pattern(Ap0, flip, permutation): 

    Ap0 stands for the data array of volume A for pattern0
    Ap1 stands for the volume A for pattern1

    The flips are applied first then the volume is transposed using 
    ext.np.transpose(X, permutation)

    Input 
    
    Ap0: ndarray ndim=3, any 3D array indeed (p0 standing for intial pattern).   
    flip: tuple (flipi, flipj, flipk) bool ex (True, True, True) is no flip (reversion of order in axis 0, 1, 2).  
    permutation: tuple of transposition ex (0, 1, 2) is identity (no transposition)

    Output

    Ap1: ndarray, same shape as Ap0, with shifts and permutation  
    
    see also orientations in nib ornt is a (n, 2) array indicating 
        in axis 0 the index of 1 to n coords in the data array, 
        in axis 1 the flip (-1) or no flip (1) to operate on the data array. 
    
    """
    
    Ap1 = Ap0.copy()
    if flip[0]: 
        Ap1[:, :, :] = Ap1[::-1, :, :]
    if flip[1]: 
        Ap1[:, :, :] = Ap1[:, ::-1, :]
    if flip[2]: 
        Ap1[:, :, :] = Ap1[:, :, ::-1]
    Ap1 = ext.np.transpose(Ap1, permutation)
    return Ap1
    
def compute_subresolution_A1(A1, aff1, nrand=100, resolution=5, fun=ext.np.mean): 
    """

    A2 = compute_subresolution_A1(A1, aff1, nrand=100, resolution=5, fun=ext.np.mean) 

    Input 
    
    A1: data array of vol1
    aff1: affine of vol1, 
    nrand=100: number of random samples for subresolution voxel
    resolution=5: resolution of the subresolution volume. resolution is in mm

    Output

    A2: ndarray with shape (IMAX, JMAX, KMAX) = (int(xmax // resolution - 1), int(ymax // resolution - 1), int(zmax // resolution - 1))
    
    
    """
    (I1, J1, K1) = A1.shape
    # print(I1, J1, K1)
    res = ext.nib.affines.apply_affine(aff1, [[0, 0, 0], [I1, J1, K1]])
    (xmax, ymax, zmax) = ext.np.abs(res[1, :] - res[0, :])
    # print(xmax, ymax, zmax)
    (IMAX, JMAX, KMAX) = (int(xmax // resolution - 1), int(ymax // resolution - 1), int(zmax // resolution - 1))
    # print(IMAX, JMAX, KMAX)
    aff1inv = ext.np.linalg.inv(aff1)
    aff2 = ext.nib.affines.rescale_affine(aff1, (I1, J1, K1), (resolution, resolution, resolution), (IMAX, JMAX, KMAX))
    # print(aff2)
    A2 = ext.np.zeros((IMAX, JMAX, KMAX))
    for i in range(IMAX): 
        print(".", end="")
        for j in range(JMAX): 
            for k in range(KMAX): 
                values = get_n_random_values_of_voxels_in_A_with_coords_in_lowres_index(
                    A1, i, j, k, I1, J1, K1, aff2, aff1inv, nrand=nrand)
                A2[i, j, k] = fun(values)
    return (A2, aff2)

def compute_subresolution_A1_aligned(
    A1, resolution1, resolution2, nrand=100, fun=ext.np.mean): 
    """
    (A2, aff2) = compute_subresolution_A1_aligned(A1, resolution1, resolution2, 
    nrand=100, fun=ext.np.mean) 

    restrictions
    the affine corresponding to A1 is consider to be 
    np.array([[s1x, 0, 0, 0], [0, s1y, 0, 0], [0, 0, s1z, 0], [0, 0, 0, 1]]) 
    the affine corresponding to A2 is consider to be 
    np.array([[s2x, 0, 0, 0], [0, s2y, 0, 0], [0, 0, s2z, 0], [0, 0, 0, 1]]) 
    We consider subsampling A1 to produce A2.
    A2 and A1 are aligned and, the center of the first voxel that is sampled is at (s2x/2, s2y/2, s2s/2). The vertex at coordinate (0,0,0) is shared with the first voxel of A1 with center at (s1x/2, s1y/2, s1z/2). 
    
    1. draw nrand sampled coordinates in each voxel of A2 and get their values in voxels of A1.
    2. compute the function on the nrand values taken in A1 and put the result in A2. 

    TODO improve code for speed or study its faisability. 
    
    
    Input

    A1: ndarray of ndim=3
    resolution1: tuple of float (resolution in mm in each dimension)
    resolution2: tuple of float (resolution in mm in each dimension)
    nrand=100: int number of samples to evaluate subsampling
    fun=ext.np.mean: function with x as input, x is an array of dim 1 containing the values to perform the function. return the statistic of the sampled values contained in x. 

    Output

    A2: nd array of ndim=3, the volume A1 at resolution2
    aff2: the affine transforms of A2, by construction contained the values of resolution2 in diagonal. 
    
    The resolution of A1 is given by
    [s1x, s1y, s1z, 1] = np.sum(A1**2, 1) ** 0.5 
    We suppose here that vectors of A1 are nearly aligned with the vectors of the standard basis.
    (I don't know what happen if this is not the case. )
    
    """
    (s1x, s1y, s1z) = resolution1
    (s2x, s2y, s2z) = resolution2
    (s1xo2, s1yo2, s1zo2) = (s1x / 2, s1y / 2, s1z / 2)
    (s2xo2, s2yo2, s2zo2) = (s2x / 2, s2y / 2, s2z / 2)
    (I1, J1, K1) = A1.shape
    (I2, J2, K2) = (int((I1 * s1x) // s2x + 1), 
                    int((J1 * s1y) // s2y + 1), 
                    int((K1 * s1z) // s2z + 1)) 
    A2 = ext.np.zeros((I2, J2, K2))
    xc2 = - s1xo2 + s2xo2
    yc2 = - s1yo2 + s2yo2
    zc2 = - s1zo2 + s2zo2
    Vijk = ext.np.zeros((nrand, ))
    for i in range(I2 - 1): 
        xi = xc2 + i * s2x
        xi1 = xi - s2xo2
        for j in range(J2 - 1): 
            yj = yc2 + j * s2y
            yj1 = yj - s2yo2
            for k in range(K2 - 1): 
                zk = zc2 + k * s2z
                zk1 = zk - s2zo2
                r = ext.np.random.rand(nrand, 3)
                Xi = xi1 + r[:, 0] * s2x
                Yj = yj1 + r[:, 1] * s2y
                Zk = zk1 + r[:, 2] * s2z
                Ii = (Xi / s1x).astype(int)
                Jj = (Yj / s1y).astype(int)
                Kk = (Zk / s1z).astype(int)
                for irand in range(nrand): 
                    Vijk[irand] = A1[Ii[irand], Jj[irand], Kk[irand]]
                A2[i, j, k] = fun(Vijk)
    aff2 = ext.np.array(
        [[s2x,   0,   0, 0], 
         [  0, s2y,   0, 0], 
         [  0,   0, s2z, 0], 
         [  0,   0,   0, 1]])
    return (A2, aff2)
    
def get_B_extract_with_translation_comparable_to_pattern(B, i, j, k, I, J, K): 
    """
    Bx = get_B_extract_with_translation_comparable_to_pattern(B, i, j, k, I, J, K): 
    extract a subarray Bx, in B at position (i,j,k) with shape (I,J,K)
    
    Input

    B: ndarray, inital array.
    i,j,k: int, start of slices where extraction is performed. 
    I,J,K: int, size of slices, will be in the shape of Bx.  

    Output

    Bx: ndarray with shape (I,J,K), the extracted array from B (comparable to a pattern, for evaluation of a similarity index or other score between Bx and the pattern Ap). 
    
    """
    
    s1 = slice(i, i + I)
    s2 = slice(j, j + J)
    s3 = slice(k, k + K)
    Bx = B[s1, s2, s3]
    return Bx

def get_range_A(A1, aff1): 
    """
    coords2 = get_range_A(A1, aff1)
    
    get the coordinates of the two points in A1 at location (0,0,0) and (I1,J1,K1)
    where (I1, J1, K1) = A1.shape

    Input

    A1: ndarray ndim=3, data of a volume
    aff1: ndarray of shape (4, 4), affine transform of the volume 

    Output

    coords2: list of float (coordinates in 3D for point (0,0,0) and point (I1, J1, K1). 
    
    
    TODO check the border at (I1, J1, K1) or (I1-1, J1-1, K1-1) ?? 
    
    
    """
    (I1, J1, K1) = A1.shape
    coords2 = ext.nib.affines.apply_affine(aff1, [[0, 0, 0], [I1, J1, K1]])
    return coords2

def get_data_from_image(A1):
    """
    A2 = get_data_from_image(A1)
    
    Input 
    
    A1: ndarray (n, n) for representation in row, col orientation for image (axis 0 for coordinate top to bottom, axis 1 for coordinate left to right). 

    return 
    
    A2: ndarray (n, n) for representation of data in i, j orientation (axis 0 for coordinate left to right, and axis 1 for coordinate bottom to top)
    """
    # same as np.rot90(A, -1) but less intuitive. 
    A2 = ext.np.flipud(A1)
    A2 = ext.np.transpose(A2)
    return A2

def get_image_from_data(A1):
    """
    A2 = get_image_from_data(A1)
    
    Input 
    
    A1: ndarray (n, n) with (axis 0 for coordinate left to right, and axis 1 for coordinate bottom to top)
    
    Output 
    
    A2: ndarray (n, n) of bool, mask (True is white and False is black) in row, col orientation for image. 
    """
    # same as np.rot90(A, -1) but less intuitive. 
    A2 = ext.np.transpose(A1)
    A2 = ext.np.flipud(A2)
    return A2

def get_values_of_voxels_in_A(A1, ijk1, I1, J1, K1):
    """
    values = get_values_of_voxels_in_A(A1, ijk1, I1, J1, K1):
    
    A is a data volume ijk1 is an array of coordinates for nvoxels. 
    get the value of the voxel or return a 0 is the voxel is not available. 
    
    Input 
    
    A1: ndarray (I,J,K) A1 is a volume
    ijk1: ndarray (nvoxels, 3) of int coordinates, 
        contains nvalues in rows with 
        ijk1[;, 0] pointing to indices in i, 
        col1 is index in j, col 2 is index in k etc. 
    (I1, J1, K1): shape of A1

    Output

    values: ndarray of shape (nvoxels, )
    
    """
    nvalues = ijk1.shape[0]
    values = ext.np.zeros((nvalues, ))
    for i in range(nvalues): 
        iA1 = ijk1[i, 0]
        jA1 = ijk1[i, 1]
        kA1 = ijk1[i, 2] 
        cond1 = (iA1 < 0) | (jA1 < 0) | (kA1 < 0)
        cond2 = (iA1 >= I1) | (jA1 >= J1) | (kA1 >= K1)
        if cond1 | cond2: 
            values[i] = 0
        else: 
            values[i] = A1[iA1, jA1, kA1].copy()
    return values 

def get_n_random_values_of_voxels_in_A_with_coords_in_lowres_index(
    A1, i2, j2, k2, I1, J1, K1, aff2, aff1inv, nrand=100): 
    """
    values = get_n_random_values_of_voxels_in_A_with_coords_in_lowres_index(A1, i2, j2, k2, I1, J1, K1, aff2, aff1inv, nrand=100)

    The sampling is done in the voxel centered at (i2, j2, k2)
    
    
    Input 

    A1: ndarray ndim=3, volume in high res where to sample
    (i2, j2, k2): int, indices of the voxel in low res
    (I1, J1, K1): int, number of elements in each axis of A1
    aff2: ndarray shape (4, 4), affine of volume A2 in lowres
    aff1inv: the inverse affine of A1, given coordinates in space, return the indices in A1. 
    nrand=100: int, number of samples 

    Output

    values: ndarray (nrand,) values sampled on A1 from voxel at (i2,j2,k2) with aff2 (low-res). 

    """
    # print(A1.shape)
    # print(aff2)
    # print(aff1inv)
    # The dimension between voxel is obtained from the voxel at (i2+1, j2+1, k2+1). 
    coords2 = ext.nib.affines.apply_affine(aff2, [[i2, j2, k2], [i2 + 1, j2 + 1, k2 + 1]])
    # print(coords1)
    r = ext.np.random.rand(nrand, 3) - 0.5 # the coord1[0, :] is the middle of the first voxel
    coords1 = ext.np.zeros((nrand, 3))
    for i in range(nrand): 
        for j in range(3): 
            coords1[i, j] = coords2[0, j] + (coords2[1, j] - coords2[0, j]) * r[i, j]
    ijk1 = ext.nib.affines.apply_affine(aff1inv, coords1).astype(int)
    values = get_values_of_voxels_in_A(A1, ijk1, I1, J1, K1)
    return values 


def remove_mask_on_A_along_axis_from_to(mask1, A1, axis, i1, i2): 
    """
    apply_mask_on_A_along_axis_from_to(mask1, A1, axis, i1, i2)

    inplace applying return A1
    
    """
    if axis == 2: 
        for k in range(i1, i2): 
            A1[:, :, k] -= mask1
    elif axis == 1:
        for j in range(i1, i2): 
            A1[:, j, :] -= mask1
    elif axis == 0:
        for k in range(i1, i2): 
            A1[k, :, :] -= mask1
    # ext.np.clip(A1, 0, 1, A1)
    return A1

def add_mask_on_A_along_axis_from_to(mask1, A1, axis, i1, i2): 
    """
    A1 = apply_mask_on_A_along_axis_from_to(mask1, A1, axis, i1, i2)
    add the value in mask on A1, at location i1 to i2 (range(i1, i2)) along axis.
    even if A1 is returned the substitution is done inplace. 

    input

    mask1: ndarray with ndim=2, indeed an image (a slice) with mask value 0 or 1. 
    A1: ndarray with shape (I,J,K)
    axis: int, the axis along which the image mask is applied (0, 1, or 2)
    i1: index of the beginning of addition of the mask
    i2: index of the end of addition (excluded). 

    output

    A1: ndarray with shape (I,J,K), A1 processed. 
    
    """
    if axis == 2: 
        for k in range(i1, i2): 
            A1[:, :, k] += mask1
    elif axis == 1:
        for j in range(i1, i2): 
            A1[:, j, :] += mask1
    elif axis == 0:
        for k in range(i1, i2): 
            A1[k, :, :] += mask1
    # ext.np.clip(A1, 0, 1, A1)
    return A1

def plot_slices(A1, mode="axial", step=1, aspect=1, decimate=1, **params): 
    """
    
    fig = plot_slices(A1, mode="axial", step=1, aspect=1, decimate=1, **params): 
    
    
    A1 is supposed to be in RAS order. (Nifti convention)
    slices are taken along axis defined by mode with step in integer values
    mode: "sagittal", "coronal", "axial"
        axis are 0, 1, 2
    LR, PA, IS
    
    input
    
    A1: ndarray with shape (I1, J1, K1)
    mode="axial": string, indicate the kind of slices to plot, select the axis alonf which to operate {"sagittal", "coronal", "axial"}.  
    step=1: int, step between slices
    aspect=1: float, determine the size of the figure (figwidth)
    decimate=1: int, if the initial figure is to big, can decimate the initial image.  
    **params: to be passed to ax.matshow (see matplotlib)

    output

    fig: matplotlib.figure.Figure
    
    """
    (I1, J1, K1) = A1.shape
    if mode == "sagittal": 
        (R, C, KMAX) = (K1, J1, I1)
    elif mode == "coronal":
        (R, C, KMAX) = (K1, I1, J1)
    elif mode == "axial":
        (R, C, KMAX) = (J1, I1, K1)
    nslices = KMAX // step
    # print("% -> " + str(KMAX % step))
    nrows = R // decimate
    ncols = C * (1 + nslices) // decimate
    # print(C)
    size = nrows * ncols
    if size >= 1000000: 
        print("! image is too big")
        return None
    B = ext.np.zeros((nrows, ncols))
    # print(I1, J1, K1, nslices)
    indices = list(range(0, KMAX, step))
    # print(B.shape, indices)
    d = decimate
    for (i, z) in enumerate(indices): 
        k1 = C * i // d
        k2 = k1 + C // d
        # print(k1, k2)
        if mode == "sagittal": 
            A2 = get_image_from_data(A1[z, ::d, ::d])
        elif mode == "coronal": 
            A2 = get_image_from_data(A1[::d, z, ::d])
        elif mode == "axial": 
            A2 = get_image_from_data(A1[::d, ::d, z])
        # print(k1, k2, A2.shape)
        B[:, k1: k2] = A2[: nrows, : k2 - k1]
    # print(k2)
    (fig, ax) = ext.plt.subplots(1, 1, figsize=(aspect * nslices, aspect))
    ax.matshow(B, **params)
    ax.axis("off")
    return fig

def compute_template_phantom_mri_highres(): 
    """
    
    (A1, aff1) = compute_template_phantom_mri_highres()
    
    
    In the mri the water raise signal in the voxel and the material does not.
    The air does not raise signal. 
    The template is then obtained by filling the tube with water and extracting the matter made by the phantom. 
    the tube is not loaded since it does not raise signal. 

    Input

    Output

    A1: ndarray with shape (I1, J1, K1), the data in the volume. 
    aff1: ndarray with shape (4, 4), the affine transform of the volume. 
    
    """
    dic1 = {}
    names = {"water", "maintainer1", "maintainer2", 
             "LRISbase", "LRIS", 
             "PAISbase", "PAIS", 
             "LRPAbase", "LRPA"}
    for name in names: 
        dic1[name] = ext.iio.imread(f"./img/{name}.pbm")
    params = [
        ("water", 1, 0, 600), 
        ("maintainer1", 1, 100, 450), 
        ("maintainer2", 1, 100, 450), 
        ("LRISbase", 1, 450, 500), 
        ("LRIS", 1, 150, 450), 
        ("LRPAbase", 2, 100, 150), 
        ("LRPA", 2, 150, 450), 
        ("PAISbase", 0, 100, 150), 
        ("PAIS", 0, 150, 450)]
    A1 = ext.np.zeros((600, 600, 600), int)
    aff1 = ext.np.array(
    [[0.1,   0,   0, 0], 
     [  0, 0.1,   0, 0], 
     [  0,   0, 0.1, 0], 
     [  0,   0,   0, 1]])
    # add water 
    mask1 = get_data_from_image(dic1[params[0][0]])
    A1 = add_mask_on_A_along_axis_from_to(mask1, A1, *params[0][1:])
    # remove others
    for param in params[1:]: 
        mask1 = get_data_from_image(dic1[param[0]])
        A1 = remove_mask_on_A_along_axis_from_to(mask1, A1, *param[1:]) 
    return (A1, aff1)

def compute_template_ct_highres(): 
    """
    
    (A1, aff1) = compute_template_phantom_mri_highres()
    
    
    In the ct the signal correspond to the mask. 
    The air does not raise signal. 
    The template is then obtained by just adding masks. 
    the tube is loaded. 

    Input

    Output

    A1: ndarray with shape (I1, J1, K1), the data in the volume. 
    aff1: ndarray with shape (4, 4), the affine transform of the volume. 
    
    """
    
    dic1 = {}
    names = {"tube", "maintainer1", "maintainer2", 
             "LRISbase", "LRIS", 
             "PAISbase", "PAIS", 
             "LRPAbase", "LRPA"}
    for name in names: 
        dic1[name] = ext.iio.imread(f"./img/{name}.pbm")
    params = [
        ("tube", 1, 0, 600), 
        ("maintainer1", 1, 100, 450), 
        ("maintainer2", 1, 100, 450), 
        ("LRISbase", 1, 450, 500), 
        ("LRIS", 1, 150, 450), 
        ("LRPAbase", 2, 100, 150), 
        ("LRPA", 2, 150, 450), 
        ("PAISbase", 0, 100, 150), 
        ("PAIS", 0, 150, 450)]
    A1 = ext.np.zeros((600, 600, 600), int)
    aff1 = ext.np.array(
    [[0.1,   0,   0, 0], 
     [  0, 0.1,   0, 0], 
     [  0,   0, 0.1, 0], 
     [  0,   0,   0, 1]])
    for param in params: 
        mask1 = get_data_from_image(dic1[param[0]])
        A1 = add_mask_on_A_along_axis_from_to(mask1, A1, *param[1:]) 
    return (A1, aff1)

def compute_affine_to_fpmat(aff1): 
    """
    A2 = compute_affine_to_fpmat(aff1)

    Compute the fpmat, the flips and permutations matrix, associated to the affine transform aff1

    Example 

    >>> aff1 = np.array([
    >>>     [0, -0.5, 0, 21], 
    >>>     [-0.25, 0, 0, 9.5], 
    >>>     [0, 0, 0.5, -24.75], 
    >>>     [0, 0, 0, 1]])
    >>> A2 = tools.compute_affine_to_fpmat(aff1)
    >>> print(A2)
    [[ 0. -1.  0.]
     [-1.  0.  0.]
     [ 0.  0.  1.]]
    
    Input 

    aff1: ndarray with shape (4, 4), the affine matrix of a volume
    
    Output

    A2: ndarray with shape (3, 3), the flips and permutation matrix with (-1) if there is a flip (1) if there is no flip, and 0 otherwise. $A2(i, j) \\neq 0$ if there is a permutation with j = p(i), 0 otherwise.   
    
    we suppose axes are approximately aligned to canonical bases
    (rotation angles are small). 
    fpmat is the matrix containing flips and permutations.
    """
    A1 = aff1[:3, :3]
    A2 = ext.np.zeros_like(A1)
    for i in range(3): 
        j = ext.np.argmax(ext.np.abs(A1[i, :]))
        sj = ext.np.sign(A1[i, j])
        A2[i, j] = sj
    return A2             

def compute_change_of_base_matrix(Xold, Xnew):
    """
    P = compute_change_of_base_matrix(Xold, Xnew):
    
    Compute change-of-base matrix from old to new base
    Xold = P @ Xnew
    
    (not really necessary in this study)

    Example

    Input

    Xold: ndarray with shape (3, 3), old basis
    Xnew: ndarray with shape (3, 3), new basis

    Output

    P: ndarray with shape (3, 3)
    
    """
    Xnewinv = ext.np.linalg.inv(Xnew)
    P = Xold @ Xnewinv
    return P 

def compute_matrix_to_affine(A1):
    """
    aff1 = compute_matrix_to_affine(A1):

    generate the affine transform from the transformation matrix (3, 3)
    A1 -> [[A1, 0], [0, 1]]

    Example 

    >>> A1 = np.arange(9).reshape((3, 3))
    >>> compute_matrix_to_affine(A1)
    array([[0., 1., 2., 0.],
           [3., 4., 5., 0.],
           [6., 7., 8., 0.],
           [0., 0., 0., 1.]])
       
    Input
    
    A1: ndarray (3, 3)

    Output
    
    aff1: ndarray (4, 4)
    
    """
    aff1 = ext.np.zeros((4, 4)) 
    aff1[:3, :3] = A1.copy()
    aff1[3, 3] = 1
    return aff1

def compute_affine_from_fp(flips, permutations): 
    """
    A2 = compute_affine_from_fp(flips, permutations): 
    
    compute affine transform associated to flips and permutations (fp)

    Example
    
    >>> tools.compute_affine_from_fp((True, False, True), (2, 1, 0))
    array([[ 0.,  0., -1.,  0.],
           [ 0.,  1.,  0.,  0.],
           [-1.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  1.]])
    
    Input
    
    flips: tuple of bool with 3 elements, (flip0, flip1, flip2) flips on axis (0, 1, 2) (True if flip, False if no flip)
    permutations: tuple of permtuations, (perm0, perm1, perm2), axis to permutate. 

    Output

    A2: ndarray with shape (4, 4)
    
    """
    A1 = ext.np.eye((4))
    A2 = ext.np.zeros((4, 4))
    for i in range(3): 
        if flips[i]: 
            # print(flips[i], A1[i, i])
            A1[i, i] = -1
        else: 
            # print(flips[i], A1[i, i])
            pass
    for i in range(3):
        j = permutations[i]
        A2[:, j] = A1[:, i].copy()
    A2[:, 3] = A1[:, 3].copy()
    return A2

def compute_transformed_affine(aff1, flips, permutations): 
    """
    aff4 = compute_transformed_affine(aff1, flips, permutations) 

    compute the new affine transform by applying th flips and permutations to aff1

    Example

    >>> aff1 = np.eye(4) 
    >>> aff1[:3, :3] *= 0.25
    >>> tools.compute_transformed_affine(aff1, (True, False, True), (2, 1, 0)) 
    array([[ 0.  ,  0.  , -0.25,  0.  ],
           [ 0.  ,  0.25,  0.  ,  0.  ],
           [-0.25,  0.  ,  0.  ,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  1.  ]])
       
    Input

    aff1: ndarray with shape (4, 4), initial affine transform
    flips: tuple of bool with 3 elements, (flip0, flip1, flip2) flips on axis (0, 1, 2) (True if flip, False if no flip)
    permutations: tuple of permtuations, (perm0, perm1, perm2), axis to permute. 

    Output

    aff4: ndarray with shape (4, 4)
    
    """
    A1 = compute_affine_to_fpmat(aff1)
    aff2 = compute_affine_from_fp(flips, permutations)
    A2 = compute_affine_to_fpmat(aff2)
    # Q = compute_change_of_base_matrix(A1, A2)
    # P = ext.np.linalg.inv(Q) # same as compute_change_of_base(A2, A1) indeed. 
    A1inv = ext.np.linalg.inv(A1)
    P = A2 @ A1inv
    aff3 = compute_matrix_to_affine(P)
    aff4 = aff3 @ aff1
    return aff4
    

def save_data_affine_to_niigz(A1, aff1, fn1): 
    """
    save_data_affine_to_niigz(A1, aff1, fn1)

    save volume in NIfTI format with extension nii.gz 

    Input

    A1: ndarray with shape (3, 3), data volume
    aff1: ndarray with shape (4, 4), affine transform of the volume
    fn1: string, filename for saving volume
    
    Output
    
    """
    
    vol1 = ext.nib.volume(A1.astype(float, aff1))
    vol1.save(fn1)
    

def compute_patterns(A0): 
    """
    patterns = compute_patterns(A0) 

    Given the data volume A0, compute all the 3! * 2³ = 6 * 8 = 48 flips and permutations configurations.

    Example 

    >>> A0 = np.arange(24).reshape(3, 2, 4)
    >>> patterns = tools.compute_patterns(A0)
    >>> patterns["False_False_False_(0, 1, 2)"]
    array([[[ 0,  1,  2,  3],
            [ 4,  5,  6,  7]],
    
           [[ 8,  9, 10, 11],
            [12, 13, 14, 15]],
    
           [[16, 17, 18, 19],
            [20, 21, 22, 23]]])
    
    >>> patterns["True_False_True_(2, 1, 0)"]
    array([[[19, 11,  3],
            [23, 15,  7]],

           [[18, 10,  2],
            [22, 14,  6]],

           [[17,  9,  1],
            [21, 13,  5]],

           [[16,  8,  0],
            [20, 12,  4]]])
    
    Input
    
    A0: ndarray with shape (I, J, K), The original pattern

    Output 

    patterns: dictionary containing the 48 patterns, 
        keys of the pattern is "s1_s2_s3_s4" 
        with 
            s1, s2, s3 in {"True", "False"}
            s4 is the permutation (p0, p1, p2) with p0, p1, p2 in {0, 1, 2}
    
    """
    Ap0 = A0.copy()
    patterns = {}
    set_permutation = {
        (0, 1, 2), 
        (0, 2, 1), 
        (1, 0, 2), 
        (1, 2, 0), 
        (2, 0, 1), 
        (2, 1, 0)}
    for ri in {True, False}: 
        for rj in {True, False}: 
            for rk in {True, False}: 
                for perm in set_permutation: 
                    s1 = str(ri)
                    s2 = str(rj)
                    s3 = str(rk)
                    s4 = str(perm)
                    k = f"{s1}_{s2}_{s3}_{s4}"
                    patterns[k] = compute_new_pattern(Ap0, (ri, rj, rk), perm)
    return patterns

def load_template(fn1): 
    """
    (A1, aff1) = load_template(fn1)

    Input

    fn1: string, filename of the template

    Output

    A1: ndarray with shape (I1, J1, K1), data of the template
    aff1: ndarray with shape (4, 4), affine transform of the template 
    
    """
    vol1 = ext.nib.load(fn1)
    A1 = vol1.get_fdata()
    aff1 = vol1.affine
    return (A1, aff1)


