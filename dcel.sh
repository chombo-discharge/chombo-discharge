#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_vertex.H\"/\#include <CD_DcelVertex.H>/g' $i
    sed -i 's/\#include <dcel_vertex.H>/\#include <CD_DcelVertex.H>/g' $i
    sed -i 's/\#include \"dcel_vertexI.H\"/\#include <CD_DcelVertexImplem.H>/g' $i
    sed -i 's/\#include <dcel_vertexI.H>/\#include <CD_DcelVertexImplem.H>/g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_edge.H\"/\#include <CD_DcelEdge.H>/g' $i
    sed -i 's/\#include <dcel_edge.H>/\#include <CD_DcelEdge.H>/g' $i
    sed -i 's/\#include \"dcel_edgeI.H\"/\#include <CD_DcelEdgeImplem.H>/g' $i
    sed -i 's/\#include <dcel_edgeI.H>/\#include <CD_DcelEdgeImplem.H>/g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_face.H\"/\#include <CD_DcelFace.H>/g' $i
    sed -i 's/\#include <dcel_face.H>/\#include <CD_DcelFace.H>/g' $i
    sed -i 's/\#include \"dcel_faceI.H\"/\#include <CD_DcelFaceImplem.H>/g' $i
    sed -i 's/\#include <dcel_faceI.H>/\#include <CD_DcelFaceImplem.H>/g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_mesh.H\"/\#include <CD_DcelMesh.H>/g' $i
    sed -i 's/\#include <dcel_mesh.H>/\#include <CD_DcelMesh.H>/g' $i
    sed -i 's/\#include \"dcel_meshI.H\"/\#include <CD_DcelMeshImplem.H>/g' $i
    sed -i 's/\#include <dcel_meshI.H>/\#include <CD_DcelMeshImplem.H>/g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_iterator.H\"/\#include <CD_DcelIterator.H>/g' $i
    sed -i 's/\#include <dcel_iterator.H>/\#include <CD_DcelIterator.H>/g' $i
    sed -i 's/\#include \"dcel_iteratorI.H\"/\#include <CD_DcelIteratorImplem.H>/g' $i
    sed -i 's/\#include <dcel_iteratorI.H>/\#include <CD_DcelIteratorImplem.H>/g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_poly.H\"/\#include <CD_DcelPoly.H>/g' $i
    sed -i 's/\#include <dcel_poly.H>/\#include <CD_DcelPoly.H>/g' $i
    sed -i 's/\#include \"dcel_polyI.H\"/\#include <CD_DcelPolyImplem.H>/g' $i
    sed -i 's/\#include <dcel_polyI.H>/\#include <CD_DcelPolyImplem.H>/g' $i
done

#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"dcel_parser.H\"/\#include <CD_DcelParser.H>/g' $i
    sed -i 's/\#include <dcel_parser.H>/\#include <CD_DcelParser.H>/g' $i
    sed -i 's/\#include \"dcel_parserI.H\"/\#include <CD_DcelParserImplem.H>/g' $i
    sed -i 's/\#include <dcel_parserI.H>/\#include <CD_DcelParserImplem.H>/g' $i
done


# Move files
mv Source/Geometry/dcel_vertex.H    Source/Geometry/CD_DcelVertex.H
mv Source/Geometry/dcel_vertexI.H   Source/Geometry/CD_DcelVertexImplem.H

mv Source/Geometry/dcel_edge.H    Source/Geometry/CD_DcelEdge.H
mv Source/Geometry/dcel_edgeI.H   Source/Geometry/CD_DcelEdgeImplem.H

mv Source/Geometry/dcel_face.H    Source/Geometry/CD_DcelFace.H
mv Source/Geometry/dcel_faceI.H   Source/Geometry/CD_DcelFaceImplem.H

mv Source/Geometry/dcel_mesh.H    Source/Geometry/CD_DcelMesh.H
mv Source/Geometry/dcel_meshI.H   Source/Geometry/CD_DcelMeshImplem.H

mv Source/Geometry/dcel_iterator.H    Source/Geometry/CD_DcelIterator.H
mv Source/Geometry/dcel_iteratorI.H   Source/Geometry/CD_DcelIteratorImplem.H

mv Source/Geometry/dcel_poly.H    Source/Geometry/CD_DcelPoly.H
mv Source/Geometry/dcel_polyI.H   Source/Geometry/CD_DcelPolyImplem.H

mv Source/Geometry/dcel_parser.H    Source/Geometry/CD_DcelParser.H
mv Source/Geometry/dcel_parserI.H   Source/Geometry/CD_DcelParserImplem.H

# for i in `find . -type f \( -iname \dcel*.H -o -iname \dcel*.cpp -o -iname \CD_Dcel*.H -o -iname \CD_Dcel*.cpp  \)`; do
#     sed -i 's/vertex/Vertex/g' $i
#     sed -i 's/edge/Edge/g' $i
#     sed -i 's/face/Face/g' $i
#     sed -i 's/mesh/Mesh/g' $i
#     sed -i 's/parser/Parser/g' $i
#     sed -i 's/iterator/Iterator/g' $i
# done
