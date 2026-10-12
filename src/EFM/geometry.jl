# Generic geometry is implemented in Jutul's cut-cell module.
import Jutul.CutCellMeshes: cross2, polygon_area_centroid, orient_polygon
import Jutul.CutCellMeshes: validate_polygon, convex_polygon, clean_polygon
import Jutul.CutCellMeshes: convex_parts, polygon_winding_parts, clip_halfplane, clip_polygon
import Jutul.CutCellMeshes: point_index!, triangle_plane_segment, cell_section
import Jutul.CutCellMeshes: merge_intervals, line_intervals, overlap_intervals
import Jutul.CutCellMeshes: polygon_average_distance, tetra_volume, tetra_average_distance
import Jutul.CutCellMeshes: cell_average_distance, project, unproject, coincident_faces, segment_overlap
import Jutul.CutCellMeshes: polygon_parts_area_centroid as patch_geometry
import Jutul.CutCellMeshes: clip_surface_segments as clipped_face_segments
import Jutul.CutCellMeshes: connected_polygon_parts as connected_patch_parts
import Jutul.CutCellMeshes: segment_union_measure as interface_measure
import Jutul.CutCellMeshes: polygon_prism as fracture_prism
import Jutul.CutCellMeshes: component_root as fracture_component_root
