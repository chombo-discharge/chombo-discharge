.. warning::

   The complex-geometry support in ``chombo-discharge`` is restricted to watertight, orientable surfaces.
   Success will fluctuate for other types of input.
   
.. _Chap:ComplexGeometries:
   
Complex geometries
------------------

Complex geometries from polygon surfaces (e.g., triangulations) are supported in ``chombo-discharge`` through half-edge mesh structures.
Since these are much more involved than analytic functions, the implementation concepts are provided in detail here.

.. note::
   
   The complex geometry support was written so that it has no external dependencies, and an extraction of the source files will therefore also compile outside of ``chombo-discharge``.

.. _Chap:DCEL:

DCEL mesh
_________

``chombo-discharge`` uses a doubly-connected edge list (DCEL) for storing polygon meshes. 
The central structure in a DCEL mesh is a half-edge, and this mesh structure is therefore also known as a "half-edge" mesh.
Half-edges circulate the inside of a polygon face, and contain references (pointers) to the *next* half-edge as well as the *pair* half-edge.
They also store a reference to the vertex where they started. 
The next half-edge is simply the next half-edge around the inside of the polygon face, whereas the pair half-edge is the half edge running in the opposite direction (which thus circulates a connected polygon face). 
Note that a polygon face is always closed, i.e. there is a finite number of half-edges that circulate the polygon.

.. figure:: /_static/figures/DCEL.png
   :width: 480px
   :align: center

   DCEL mesh structure. Each half-edge stores references to previous/next half-edges, the pair edge, and the starting vertex.
   Vertices store a coordinate as well as a reference to one of the outgoing half-edges.

In ``chombo-discharge``, the polygon face class stores a reference to one of the half-edges circulating the polygon.
Since this half-edge can be used to build a list of half-edges that circulate the polygon, as well as a list of vertices that make up the polygon face, one can easily compute the face area, centroid, and normal vector.
For efficiency reasons, these quantities are stored for each polygon face. 
Note that polygon faces are not restricted to triangles.

.. note::

   The DCEL implementation in ``chombo-discharge`` is found in :file:`$DISCHARGE_HOME/Source/Geometry` and consists of all files name :file:`CD_Dcel*`.

Vec2T and Vec3T
^^^^^^^^^^^^^^^

The polygon surface geometry benefits from decoupling from ``Chombo`` dimensionality and precision.
The reason for this is:

#. DCEL grids consume memory, but typically do not require double precision. 
#. DCEL grids are always in 3D, but ``chombo-discharge`` can also run in 2D.

``chombo-discharge`` uses always-2D and always-3D vector types with template precision for this purpose, and these are precision template as ``Vec2T<T>`` and ``Vec3T<T>``. 
This is done in order to reduce memory, as well as support use of polygon surfaces in 2D simulations. 
In this case one can obtain a 2D slide of the 3D surface by setting the :math:`z` value to some specified value. 

.. _Chap:DCELDistance:

Signed distance
^^^^^^^^^^^^^^^

To compute the signed distance from a point :math:`\mathbf{x}` to a DCEL mesh, one must be able to compute the signed distance to vertices, edges, or polygon faces. 
Which feature is closest depends on the projection of the point :math:`\mathbf{x}` to the plane of the polygon face, see :numref:`Fig:PolygonProjection`. 

.. _Fig:PolygonProjection:
.. figure:: /_static/figures/PolygonProjection.png
   :width: 240px
   :align: center

   Possible closest-feature cases after projecting a point :math:`\mathbf{x}` to the plane of a polygon face. 

The signed distance function is therefore checked as follows:

#. **Polygon face**.
   a
   When computing the distance from a point :math:`\mathbf{x}` to the polygon face we first determine if the projection of :math:`\mathbf{x}` to the face's plane lies inside or outside the face.
   This is more involved that one might think, and in ``chombo-discharge`` this is done by first computing the two-dimensional projection of the polygon face, ignoring one of the coordinates.
   Next, we determine, using 2D algorithms, if the projected point lies inside the embedded 2D representation of the polygon face. 
   Various algorithms for this are available, such as computing the winding number, the crossing number, or the subtended angle between the point and the 2D polygon.

   .. note::
   
      ``chombo-discharge`` will use the crossing number algorithm by default.
      
   If the point projects to the inside of the face, the signed distance is just :math:`d = \mathbf{n}_f\cdot\left(\mathbf{x} - \mathbf{x}_f\right)` where :math:`\mathbf{n}_f` is the face normal and :math:`x_f` is a point on the face plane (e.g., a vertex).
   Otherwise, the closest feature is an edge or a vertex.
   
#. **Edge**.
   
   When computing the signed distance to an edge, the edge is parametrized as :math:`\mathbf{e}(t) = \mathbf{x}_0 + \left(\mathbf{x}_1 - \mathbf{x}_0\right)t`, where :math:`\mathbf{x}_0` and :math:`\mathbf{x}_1` are the starting and ending vertex coordinates.
   The point :math:`\mathbf{x}` is projected to this line, and if the projection yields :math:`t^\prime \in [0,1]` then the edge is the closest point.
   In that case the signed distance is the projected distance and the sign is given by the sign of :math:`\mathbf{n}_e\cdot\left(\mathbf{x} - \mathbf{x}_0\right)` where :math:`\mathbf{n}_e` is the pseudonormal vector of the edge. 
   Otherwise, the closest point is one of the vertices.
   
#. **Vertex**.

   If the closest point is a vertex then the signed distance is simply :math:`\mathbf{n}_v\cdot\left(\mathbf{x}-\mathbf{x}_v\right)` where :math:`\mathbf{n}_v` is the vertex pseudonormal and :math:`\mathbf{x}_v` is the vertex position.

The above signed distance computations require the use of normal vectors for faces and vertices.
We use the pseudonormal from :cite:`1407857` where an edge normal is given by

.. math::

   \mathbf{n}_{e} = \frac{1}{2}\left(\mathbf{n}_{f} + \mathbf{n}_{f^\prime}\right).

where :math:`f` and :math:`f^\prime` are the two faces connecting the edge.
The vertex pseudonormal are given by

.. math::

  \mathbf{n}_{v} = \frac{\sum_i\alpha_i\mathbf{n}_{f_i}}{\left|\sum_i\alpha_i\right|},

where the sum runs over all faces which share :math:`v` as a vertex, and where :math:`\alpha_i` is the subtended angle of the face :math:`f_i`, see :numref:`Fig:Pseudonormal`. 

.. _Fig:Pseudonormal:
.. figure:: /_static/figures/Pseudonormal.png
   :width: 240px
   :align: center

   Edge and vertex pseudonormals.

When computing the signed distance to a DCEL mesh, one would iterate through the faces of the mesh.
In ``chombo-discharge`` this would look like the following:

.. code-block:: c++

   Vec3T<double> x;
   Dcel::MeshT<double> mesh;

   double shortestDistance = std::numeric_limits<Real>::infinity();

   const std::vector<std::shared_ptr<FaceT<double> > >& faces = mesh->getFaces();

   for (const auto& f : faces){
      const double curDistance = f->signedDistance(x);

      if(std::abs(curDistance) < std::abs(shortestDistance)){
         shortestDistance = curDistance;
      }
   }

Although valid code, iterating through faces this way can be extremely time-consuming.
A DCEL mesh can contain tens of thousands of faces, and for each face we need to determine if the input point projects to the inside or outside of the polygon face, determine if the point projects to an edge, or compute the distance to a vertex.
This is to be compared to e.g. the cost of computing the signed distance for an analytic function like :math:`d(\mathbf{x}) = R - \left|\mathbf{x}\right|`. 
For this reason the DCEL functionality is almost always used together with accelerating data structures that speeds up this process by orders of magnitude (see :ref:`Chap:BVH`).

.. _Chap:DCELParser:

File parsers
____________

The DCEL functionality in ``chombo-discharge`` can support arbitrary-sided polygon faces, and mixed-polygon surface tesselations.
Generating such grids is done by reading an input file which builds the DCEL mesh.
The functionality for this is placed in :file:`$DISCHARGE_HOME/Source/Geometry/CD_DcelParser.H`.
Users are themselves responsible for ensuring that the mesh is provided as a watertight, orientable surface.
``chombo-discharge`` makes no attempt at repairing input grids, although it will issue warnings if the input mesh has holes or incomplete faces. 

.. note::

   Currently, DCEL parsing functionality is limited to reading non-binary PLY files. 
   Contributions of file parsers for other types of files are highly welcome. 

.. _Chap:BVH:

Bounding volume hierarchy
_________________________

Bounding Volume Hierarchies (BVHs) are comparatively simple data structures that can accelerate closest-point searches by orders of magnitude.
BVHs are tree structures where the regular nodes are bounding volumes that enclose all geometric primitives (e.g. polygon faces) further down in the hierarchy.
The *leaf nodes* contain a subset of the primitives, as well as a bounding volume for enclosing them, see :numref:`Fig:TrianglesBVH`.

.. note::

   Currently, spheres and axis-aligned bounding boxes (AABBs) are supported as bounding volumes.
   
To accelerate the signed distance function queries, we use a top-down constructed BVH with bounding volumes for enclosing sets of polygon faces. 
The top-down construction begins by placing all the primitives in the root node of the tree, and then partitioning the primitives into two subsets containing approximately half of the primitives each.
An example of this is shown in :numref:`Fig:TrianglesBVH` where the primitives in a parent node, whose bounding volume is grown for visual clarity, are split into left and right leaf nodes.

.. _Fig:TrianglesBVH:
.. figure:: /_static/figures/TrianglesBVH.png
   :width: 480px
   :align: center

   Example of BVH partitioning for enclosing triangles. The regular node :math:`P` contains two leaf nodes :math:`L` and :math:`R` which contain the primitives (triangles).

Bounding volumes
^^^^^^^^^^^^^^^^

Currently, two types of bounding volumes are supported: Spheres and axis-aligned bounding boxes.
These two are implemented in :file:`$DISCHARGE_HOME/Source/Geometry/CD_BoundingVolumes.H` and they contain the necessary functionality for operating with our BVH implemetation.
Note that the bounding volume is a template parameter to the BVH implementation, and that users can supply need bounding volumes if they have them. 

Tree pruning
^^^^^^^^^^^^

For :math:`N` primitives, the complexity of a direct search is :math:`\mathcal{O}\left(N\right)` whereas for a bounding volume hierarchy it can be made :math:`\mathcal{O}\left(\log N\right)`.
The reduction in complexity is enabled by pruning branches in the bounding volume tree.

The search for the closest primitive :math:`P_c` begins at the root node of the tree.
Since the primitives are stored in the leaf nodes, the search starts by selecting one of the sub-trees to investigate.
The implementation in ``chombo-discharge`` supports two algorithms here (typically, the first approach is faster since it statistically tends to prune more branches than the second)

#. First search the sub-tree whose bounding volume is closest to :math:`\mathbf{x}`. 
#. Always search the left tree first

.. note::

   For efficiency, ``chombo-discharge`` computes the unsigned square distance rather than the signed square distance in the pruning process.

This search process then recurses; at each node we select a new sub-tree to search through, which continues until a leaf node is reached. 
At the leaf node the distance to the primitives is computed directly.
As more sub-trees are subsequently investigated, the sub-trees can be pruned if the distance to the node's bounding volume is larger than the shortest distance found so far.
This process can be used to prune entire subtrees from the search, see :numref:`Fig:TreePruning`.

.. _Fig:TreePruning:
.. figure:: /_static/figures/TreePruning.png
   :width: 480px
   :align: center

   Example of BVH tree pruning. When descending from node :math:`P` we determine that we first investigate the left subtree (node :math:`A`) since its bounding volume is closer than the bounding volumes for the other subtree.
   Since :math:`A` is a leaf node, we find the signed distance from :math:`\mathbf{x}` to the primitives in :math:`A`.
   This distance is shorter than the distance from :math:`\mathbf{x}` to the bounding volume that encloses nodes :math:`B` and :math:`C`.
   Furthermore, the distance from :math:`\mathbf{x}` to any primitive in :math:`B` or :math:`C` is equal to or longer than the distance to the bounding volumes themselves.
   Consequently, the entire subtree containing :math:`B` and :math:`C` can be pruned. 


Implementation
^^^^^^^^^^^^^^

In ``chombo-discharge`` the BVH functionality for geometries is encapsulated by a template class ``NodeT<T, P, BV>``.
This class is found in :file:`$DISCHARGE_HOME/Source/Geometry/CD_BVH.H`, and implemented in :file:`$DISCHARGE_HOME/Source/Geometry/CD_BVHImplem.H`.
The interpretation of the template parameters are as follows:

* ``T`` Floating-point precision.
* ``P`` Primitive type. For DCEL functionality this will be the DCEL faces, i.e., ``DCEL::FaceT<T>``.
* ``BV`` Bounding volume type, either bounding spheres or AABBs. 

The reason for the template is that we wish to decouple the precision from that in ``Chombo``, be able to use different kinds of primitives (for future-proofing code), and to be able to use different bounding volumes (for performance reasons).

The BVH functionality uses feature-rich code, with a lot of functionality captured by a comparatively few lines of C++ code. 
To build the tree, one begins by creating a single leaf node containing all the primitives (i.e. the root node).
After that, calling ``NodeT<T, P, BV>::topDownSortAndPartitionPrimitives`` will use the top-down sorted primitive partitioning.
The input arguments to this routine are polymorphic lambdas that provide the opportunity to use user-supplied criteria for

#. Terminating the partitioning process
#. Deciding how to partition the sub-trees.
#. Construct a bounding volume from a set of primitives.

The function signature for ``topDownSortAndPartitionPrimitives`` is

.. code-block:: c++

    void topDownSortAndPartitionPrimitives(const StopFunction&      a_stopFunc,
					   const PartitionFunction& a_partFunc,
					   const BVConstructor&     a_bvFunc) noexcept;

where the function arguments are aliases for

.. code-block:: c++

   // I.e., a dynamic list of pointers to immutable primitives
   template <class P>
   using PrimitiveListT = std::vector<std::shared_ptr<const P> >;		

   // I.e., a function which returns true if the input node can be partitioned further
   template <class T, class P, class BV>
   using StopFunctionT = std::function<bool(const NodeT<T, P, BV>& a_node)>;

   // I.e. a function which splits a list of primitives into two new lists. These
   // are the primitives contained in the left/right child nodes in the tree. 
   template <class P>
   using PartitionFunctionT = std::function<std::pair<PrimitiveListT<P>, PrimitiveListT<P> >(const PrimitiveListT<P>& a_primitives)>;

   // I.e. a function which constructs a bounding volume from a list of primitives.
   template <class P, class BV>
   using BVConstructorT = std::function<BV(const PrimitiveListT<P>& a_primitives)>;

For further documentation, see :file:`$DISCHARGE_HOME/Source/Geometry/CD_BVH.H`.

.. note::

   Default implementations of these for DCEL functionality are available in :file:`$DISCHARGE_HOME/Source/Geometry/CD_DcelBVH.H`.

Implicit functions
__________________

There are two implicit functions that encapsulate the above functionality in ``chombo-discharge``.

#. ``DcelSdf`` (see :file:`$DISCHARGE_HOME/ImplicitFunctions/CD_DcelSdf.H`).
   This class does not use a bounding volume hierarchy, and the signed distance function looks through every face in the DCEL mesh.

#. ``BvhSdf`` (see :file:`$DISCHARGE_HOME/ImplicitFunctions/CD_BvhSdf.H`).
   This class uses a bounding volume hierarchy for enclosing DCEL faces.
   The distance function uses the BVH accelerating structure, and is therefore several orders of magnitude faster than ``DcelSdf``. 

Users that are looking to use DCEL and BVH functionality for describing objects should use these two implicit functions directly.
Note that both objects take a DCEL mesh in their constructor arguments, so the mesh must be constructed by parsing the input mesh from file (see :ref:`Chap:DCELParser`).
Various geometries have been tested using these classes (see :numref:`Fig:PLYObjects`). 

.. _Fig:PLYObjects:
.. figure:: /_static/figures/PLYObjects.png
   :width: 480px
   :align: center

   Examples of complex geometries built using ``chombo-discharge``.

Note that the signed-distance functions that are constructed from DCEL mesh structures act *precisely* as any other signed distance function.
That is, they can be translated, rotated, cut, or useed in CSG just like any other implicit function. 
